

from fenics import sqrt
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
        self.wind_speed_x = 0.
        self.wind_speed_y = 0.
        self.pressure = 0.
        self.aux = []
        self.aux_norm = 0.
        self.deltax = 0.
        self.wvel = 0.

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

        # get auxiliary point for wind drag theta calculation
        self.deltax = - self.tvy_curr / self.tvx_curr
        self.aux = [-1.0, self.deltax]
        self.aux_norm = sqrt(self.aux[0] ** 2. + self.aux[1] ** 2.)

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

        self.wind_speed_x = wind_speed_x
        self.wind_speed_y = wind_speed_y
        self.wvel = sqrt(self.wind_speed_x ** 2 + self.wind_speed_y ** 2)
        self.pressure = pressure

        return wind_speed_x, wind_speed_y, pressure

    def get_wind_drag(self, x_coord_deg, y_coord_deg):
        # determine which sector it is.
        cen = [x_coord_deg - self.lon_curr, y_coord_deg - self.lat_curr]
        cen_norm = sqrt(cen[0] ** 2 + cen[1] ** 2)
        if cen_norm != 0:
            theta = np.degrees(np.arccos(np.dot(self.aux, cen) / self.aux_norm / cen_norm))
        else:
            theta = 0.0

        # determin which sector it is.
        sector = None
        if (self.deltax >= 0 and y_coord_deg + self.deltax * x_coord_deg >= self.lat_curr + self.deltax * self.lon_curr) \
            or (self.deltax < 0 and y_coord_deg + self.deltax * x_coord_deg <= self.lat_curr + self.deltax * self.lon_curr):
            if 0. <= theta < 20.:
                sector = "left"
            elif 20. <= theta < 150.:
                sector = "right"
            else:
                sector = "rear"
        else:
            if 0. <= theta < 120.:
                sector = "left"
            else:
                sector = "rear"

        # determine the wind drag point by point.
        if 0. <= self.wvel <= 18.:
            wdrag = 0.001 * (0.75 + 40. / 600. * self.wvel)
        else:
            if sector == "left":
                if 18. < self.wvel <= 25.:
                    wdrag = 0.0018
                elif 25. < self.wvel <= 30.:
                    wdrag = 0.00054 * self.wvel - 0.0117
                elif 30. < self.wvel <= 45:
                    wdrag = -0.0035 / 15. * self.wvel + 0.0115
                elif self.wvel > 45.:
                    wdrag = 0.001
            elif sector == "right":
                if 18. < self.wvel <= 18.7:
                    wdrag = 0.001 * (0.75 + 40. / 600. * self.wvel)
                elif 18.7 < self.wvel <= 35.:
                    wdrag = 0.002
                elif 35. < self.wvel <= 45.:
                    wdrag = 0.0001 * self.wvel - 0.0015
                elif self.wvel > 45.:
                    wdrag = 0.003
            elif sector == "rear":
                if 18. < self.wvel <= 18.7:
                    wdrag = 0.001 * (0.75 + 40. / 600. * self.wvel)
                elif 18.7 < self.wvel <= 35.:
                    wdrag = 0.002
                elif 35. < self.wvel <= 45.:
                    wdrag = -0.0001 * self.wvel + 0.0055
                elif self.wvel > 45.:
                    wdrag = 0.001

        return wdrag




        





















