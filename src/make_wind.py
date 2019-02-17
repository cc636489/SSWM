

from fenics import sqrt, Mesh, FunctionSpace, VectorFunctionSpace, Expression, File, Function, project, as_vector
from netCDF4 import Dataset
from datetime import datetime, timedelta
import numpy as np


def get_distance(r_earth, deg_2_rad, ddx, ddy, y1, y2):
    distance = r_earth * 2 * np.arcsin(
        sqrt((np.sin(ddy / 2)) ** 2 + np.cos(y1 * deg_2_rad) * np.cos(y2 * deg_2_rad) * (np.sin(ddx / 2)) ** 2))
    return distance


class MakeWind:

    def __init__(self, inputs):

        self.inputs = inputs

        self.tim = []  # datetime object
        self.lat = []
        self.lon = []
        self.cpr = []
        self.spd = []
        self.rmw = []
        self.tvx = []
        self.tvy = []
        self.lat_curr = 0.
        self.lon_curr = 0.
        self.spd_curr = 0.
        self.cpr_curr = 0.
        self.rmw_curr = 0.
        self.tvx_curr = 0.
        self.tvy_curr = 0.
        self.delta_press = 0.
        self.ts = 0.
        self.bb = 0.

        self._get_wind_p1()
        self._get_wind_p2()

    def _get_wind_p1(self):
        with open(self.inputs.input_dir + self.inputs.wind_file) as f_wind:
            for this_line in f_wind:
                words = this_line.split()
                self.tim.append(datetime.strptime(words[2], '%Y%m%d%H'))
                self.lat.append(float(words[6]) / 10)
                self.lon.append(-float(words[7]) / 10)
                self.spd.append(float(words[8]) * 0.51444444)
                self.cpr.append(float(words[9]) * 100)
                self.rmw.append(float(words[20]) * 1.85200000318079 * 1000)
        return self

    def _get_wind_p2(self):
        for m in range(len(self.lat) - 1):
            avg_lat = (self.lat[m] + self.lat[m + 1]) / 2
            dx = get_distance(self.inputs.R_earth, self.inputs.DEG2RAD, (self.lon[m + 1] - self.lon[m]) *
                              self.inputs.DEG2RAD, 0.0, avg_lat, avg_lat)
            dy = get_distance(self.inputs.R_earth, self.inputs.DEG2RAD, 0.0, (self.lat[m + 1] - self.lat[m]) *
                              self.inputs.DEG2RAD, self.lat[m], self.lat[m + 1])
            u_trans = np.sign(self.lon[m + 1] - self.lon[m]) * dx / self.inputs.wind_dt
            v_trans = np.sign(self.lat[m + 1] - self.lat[m]) * dy / self.inputs.wind_dt
            if m == 0:
                self.tvx.append(u_trans)
                self.tvy.append(v_trans)
            self.tvx.append(u_trans)
            self.tvy.append(v_trans)
        return self

    def _get_wind_p_curr(self, current_time):
        curr_time_obj = self.tim[0] + timedelta(seconds=current_time)
        for n in range(len(self.lat) - 1):
            if (curr_time_obj >= self.tim[n]) and (curr_time_obj <= self.tim[n + 1]):
                ratio = float((curr_time_obj - self.tim[n]).seconds) / float((self.tim[n + 1] - self.tim[n]).seconds)
                self.lat_curr = ratio * (self.lat[n + 1] - self.lat[n]) + self.lat[n]
                self.lon_curr = ratio * (self.lon[n + 1] - self.lon[n]) + self.lon[n]
                self.spd_curr = ratio * (self.spd[n + 1] - self.spd[n]) + self.spd[n]
                self.cpr_curr = ratio * (self.cpr[n + 1] - self.cpr[n]) + self.cpr[n]
                self.rmw_curr = ratio * (self.rmw[n + 1] - self.rmw[n]) + self.rmw[n]
                self.tvx_curr = ratio * (self.tvx[n + 1] - self.tvx[n]) + self.tvx[n]
                self.tvy_curr = ratio * (self.tvy[n + 1] - self.tvy[n]) + self.tvy[n]
                # print curr_time_obj
                # print ratio
                # print self.lat_curr
        return self

    def get_prepared(self, current_time):
        # get current interpolated wind parameter.
        self._get_wind_p_curr(current_time)

        # get central pressure deficit.
        self.delta_press = self.inputs.P_background * 100 - self.cpr_curr
        if self.delta_press < 100.0:
            self.delta_press = 100.

        # get translational speed.
        self.ts = np.sqrt(self.tvx_curr ** 2.0 + self.tvy_curr ** 2.0)
        self.spd_curr = self.spd_curr - self.ts
        self.spd_curr = self.spd_curr / self.inputs.bla_dj
        self.bb = self.inputs.rho_air * np.e * self.spd_curr ** 2.0 / self.delta_press
        if self.bb < 1.0:
            self.bb = 1.0
        if self.bb > 2.5:
            self.bb = 2.5

    def get_wind_x_wind_y(self, x_coord_deg, y_coord_deg):
        # get relative distance from hurricane eye.
        dx_rad = (x_coord_deg - self.lon_curr) * self.inputs.DEG2RAD
        dy_rad = (y_coord_deg - self.lat_curr) * self.inputs.DEG2RAD
        theta_angle = np.arctan2(dy_rad, dx_rad)
        rad = get_distance(self.inputs.R_earth, self.inputs.DEG2RAD, dx_rad, dy_rad, self.lat_curr, y_coord_deg)
        coriolis = 2 * self.inputs.OMEGA * np.sin(y_coord_deg * self.inputs.DEG2RAD)

        # get surface pressure(Pa), wind_x(m/s), wind_y(m/s)
        pressure = self.cpr_curr + self.delta_press * np.exp(-(self.rmw_curr / rad) ** self.bb)
        vr = np.sqrt((self.rmw_curr / rad) ** self.bb * np.exp(1 - (self.rmw_curr / rad) ** self.bb) *
                     self.spd_curr ** 2 + rad ** 2 * coriolis ** 2 / 4.0) - rad * coriolis / 2.0

        # get trans_speed.
        trans_spd_x = abs(vr) / self.spd_curr * self.tvx_curr
        trans_spd_y = abs(vr) / self.spd_curr * self.tvy_curr
        # print "tvx: ", self.tvx_curr, "  tvy: ", self.tvy_curr

        # get rotational wind_x wind_y.
        wind_speed_x = - abs(vr) * np.sin(theta_angle)
        wind_speed_y = abs(vr) * np.cos(theta_angle)

        # get wind speed from surface to 10meter, and 3second gust speed to 10min average speed.
        wind_speed_x = wind_speed_x * self.inputs.bla_dj * self.inputs.one2ten
        wind_speed_y = wind_speed_y * self.inputs.bla_dj * self.inputs.one2ten

        # get total wind speed at this node at curr_time.
        wind_speed_x = wind_speed_x + trans_spd_x
        wind_speed_y = wind_speed_y + trans_spd_y

        return wind_speed_x, wind_speed_y, pressure


def test1():

    file_dir = '/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-github/input/'
    out_dir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/wind_ike/'
    wind_file_name = 'fort.22_ike_short_short'
    mesh_file_name = 'fort.14_small_gulf.xml'
    bath_file_name = 'bathymetry_small_gulf.nc'
    wind_output_file_name = "wind_xy_func_ike_short_short.pvd"
    days = 2.5

    this_mesh = Mesh(file_dir + mesh_file_name)
    this_scalar_space = FunctionSpace(this_mesh, "CG", 1)
    this_vector_space = VectorFunctionSpace(this_mesh, "CG", degree=1, dim=2)
    code = '''

    namespace dolfin {

    struct node {
        double x;
        double y;
        double v;
        node *next;
    };

    class MyFun : public Expression
    {
        private:
            node *root;

        public:

            MyFun(): Expression()
            {
                root = new node;
                root->next = 0;
            };

            void initial(double _x, double _y, double _v)
            {
                node *newNode = new node;
                newNode->x = _x;
                newNode->y = _y;
                newNode->v = _v;
                newNode->next = root;
                root = newNode;
                //cout << _v << endl;
            };

            void update(double _x, double _y, double _v)
            {
                node *p = findp(_x, _y);
                //cout << p->v << "  " << _v << endl;
                p->v = _v;
            };

            void eval(Array<double>& values, const Array<double>& x) const
            {
                double vv = findval(x[0], x[1]);
                values[0] = vv;
                //cout << x[0] << "  " << x[1] << "  " << vv << endl;
                //cout << "  " << endl;
            };

            node * findp(double _x, double _y) const
            {
                node *p = root;
                while (p->next != 0)
                {
                    if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                    {
                        return p;
                    }
                    else
                    {
                        p = p->next;
                    }
                }
                return 0;
            }

            double findval(double _x, double _y) const
            {
                node *p = root;
                while (p->next != 0)
                {   
                    if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                    {
                        //cout << fabs(p->x-_x) << "  " << fabs(p->y-_y) << endl;
                        double find = p->v;
                        return find;
                    }
                    else
                    {
                        p = p->next;
                    }
                }
                return 0;
            };

    };
    };

    '''

    # create a list for wind_para_x, wind_para_y.
    wind_para_x = Expression(cppcode=code, element=this_scalar_space.ufl_element())
    wind_para_y = Expression(cppcode=code, element=this_scalar_space.ufl_element())

    # read in variables.
    nc = Dataset(file_dir + bath_file_name, 'r', format='NETCDF4')
    x_coord = nc.variables['x_coord'][:].data.tolist()
    y_coord = nc.variables['y_coord'][:].data.tolist()
    x_deg = nc.variables['lon'][:].data.tolist()
    y_deg = nc.variables['lat'][:].data.tolist()
    bath_y = nc.variables['bathy'][:].data.tolist()

    # initialize a list for wind_para_x, wind_para_y.
    for st in range(len(x_coord)):
        wind_para_x.initial(x_coord[st], y_coord[st], bath_y[st])
        wind_para_y.initial(x_coord[st], y_coord[st], bath_y[st])

    # now we can update the wind speed in every time step.
    curr_time = 0
    final_time = 3600 * 24 * days
    dt = 3600

    ######################################################
    # this is a bug, because "MakeWind" has been modified.
    ######################################################
    wind = MakeWind(file_dir + wind_file_name)
    wind_x_max = []
    wind_y_max = []
    wind_xy_output_file = File(out_dir + wind_output_file_name)
    w = Function(this_vector_space, name="wind_speed")

    while curr_time <= final_time:

        wind.get_prepared(current_time=curr_time)
        wind_para_x_list = []
        wind_para_y_list = []
        for t in range(len(x_deg)):
            wind_x, wind_y, pressure = wind.get_wind_x_wind_y(x_coord_deg=x_deg[t], y_coord_deg=y_deg[t])
            wind_para_x_list.append(wind_x)
            wind_para_y_list.append(wind_y)

        wind_x_max.append(np.max(wind_para_x_list))
        wind_y_max.append(np.max(wind_para_y_list))

        for sm in range(len(x_coord)):
            wind_para_x.update(x_coord[sm], y_coord[sm], wind_para_x_list[sm])
            wind_para_y.update(x_coord[sm], y_coord[sm], wind_para_y_list[sm])

        # write to file
        temp = project(as_vector([wind_para_x, wind_para_y]), this_vector_space)
        w.assign(temp)
        wind_xy_output_file << (w, float(curr_time))

        # update current time
        print curr_time/dt
        curr_time += dt

    return 1


def main():

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test1()

    if flag:
        print blue + bold + "Test1: " + bold + green + "PASS"
    else:
        print blue + bold + "Test1: " + bold + red + "FAIL"

    return 0


if __name__ == "__main__":
    main()
