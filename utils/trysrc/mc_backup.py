# Start calculating MC( order = 0 ) here. Will use global variables. but will not change it inside for sure.
# def func_mc(sample_vs, sample_fric, sample_f_x, sample_f_y, sample_bathymetry, sample_ic_u, sample_ic_eta, sample_bc_u, sample_bc_eta):
#
#     t_mc = Constant(start_time)
#     basis_struct_mc = MakeStoBasis(dist_name, 0, sto_poly_dim, coefficient)
#
#     # Get horizontal domain setting
#     W_mc = VectorFunctionSpace(mesh, "CG", 2, dim=2)
#     P_mc = FunctionSpace(mesh, "CG", 1)
#
#     # Get stochastic coefficient
#     nu_expression_list_mc, friction_expression_list_mc, f_u_expression_list_mc, f_v_expression_list_mc = \
#         MakeStoModes(basis_struct_mc, sample_vs, sample_fric, sample_f_x, sample_f_y)
#
#     # Get the bathymetry
#     project_bottom_function_mc = MakeBathy(sample_bathymetry, basis_struct_mc, P_mc)
#
#     # Get initial condition in stochastic setting
#     u_ic_function_mc, eta_ic_function_mc = MakeInitialObjectFunc(sample_ic_u, sample_ic_eta, basis_struct_mc, W_mc, P_mc)
#
#     # Get boundary condition in stochastic setting
#     u_bc_object_list_mc, eta_bc_object_list_mc, uTimeDependent_mc, etaTimeDependent_mc, u_list_expression_mc, eta_list_expression_mc = MakeBoundaryObjectList(sample_bc_u, sample_bc_eta, basis_struct_mc, bc_file, mesh, domain, W_mc, P_mc, t_mc, project_bottom_function_mc)
#
#     ##################################################################
#     # starting initializing
#     ##################################################################
#
#     # set up solution function
#     v_mc = TestFunction(W_mc)
#     u_mc = TrialFunction(W_mc)
#     q_mc = TestFunction(P_mc)
#     eta_mc = TrialFunction(P_mc)
#
#     u00_mc = Function(W_mc)
#     u0_mc = Function(W_mc, name="u0")
#     ut_mc = Function(W_mc) # Tentative velocity
#     u1_mc = Function(W_mc, name="u")
#     eta0_mc = Function(P_mc, name="eta0")
#     eta1_mc = Function(P_mc, name="eta")
#
#     # define total water column
#     if linear_divergence:
#         H_mc = project_bottom_function_mc
#     else:
#         H_mc = project_bottom_function_mc + eta0_mc
#
#     # Assign initial condition
#     u0_mc.assign(u_ic_function_mc)
#     u00_mc.assign(u_ic_function_mc)
#     u1_mc.assign(u_ic_function_mc)
#     eta0_mc.assign(eta_ic_function_mc)
#     eta1_mc.assign(eta_ic_function_mc)
#
#     # about fu fv and nu
#     if "t" in f_u_expression_list_mc[0]:
#         f_u_expression_object_list_mc = Expression(f_u_expression_list_mc[0], element=W_mc.sub(0).ufl_element(), t=t_mc)
#     else:
#         f_u_expression_object_list_mc = Expression(f_u_expression_list_mc[0], element=W_mc.sub(0).ufl_element())
#     if "t" in f_v_expression_list_mc[0]:
#         f_v_expression_object_list_mc = Expression(f_v_expression_list_mc[0], element=W_mc.sub(1).ufl_element(), t=t_mc)
#     else:
#         f_v_expression_object_list_mc = Expression(f_v_expression_list_mc[0], element=W_mc.sub(1).ufl_element())
#     friction_expression_object_list_mc = Expression(friction_expression_list_mc[0], element=P_mc.ufl_element())
#     nu_expression_mc = Expression(nu_expression_list_mc[0], element=P_mc.ufl_element())
#
#     if include_les:
#         les_V_mc = FunctionSpace(mesh, "CG", 1)
#         les_mc = LES(les_V_mc, u0_mc, les_parameters['smagorinsky_coefficient'])
#         nu_mc = nu_expression_mc + les_mc.eddy_viscosity
#     else:
#         nu_mc = nu_expression_mc
#
#
#     # Tentative Velocity step
#     u_mean_mc = theta * u_mc + (1. - theta) * u0_mc
#     u_bash_mc = 3./2 * u0_mc - 1./2 * u00_mc
#     u_diff_mc = u_mc - u0_mc
#
#
#     F_u_tent_mc = (1/dt) * inner(v_mc, u_diff_mc) * dx + g * inner(v_mc, grad(eta0_mc)) * dx \
#         + friction_expression_object_list_mc * inner(u_mean_mc, v_mc) * dx \
#         - v_mc[0] * f_u_expression_object_list_mc * dx - v_mc[1] * f_v_expression_object_list_mc * dx
#
#     if include_viscosity:
#         F_u_tent_mc += nu_mc * inner(grad(v_mc), grad(u_mean_mc)) * dx
#
#     if include_convection:
#         F_u_tent_mc += inner(v_mc, grad(u_mean_mc)*u_bash_mc) * dx
#
#     a_u_tent_mc = lhs(F_u_tent_mc)
#     L_u_tent_mc = rhs(F_u_tent_mc)
#
#
#
#     # Pressure correction step
#     eta_diff_mc = eta_mc - eta0_mc
#     ut_mean_mc = theta * ut_mc + (1. - theta) * u0_mc
#     F_p_corr_mc = eta_diff_mc*q_mc * dx + theta**2 * g * dt**2 * H_mc * inner(grad(q_mc), grad(eta_diff_mc)) * dx \
#         + dt * q_mc * div(H_mc*ut_mean_mc) * dx
#     a_p_corr_mc = lhs(F_p_corr_mc)
#     L_p_corr_mc = rhs(F_p_corr_mc)
#
#
#
#     # Velocity correction step
#     eta_diff_mc = eta1_mc - eta0_mc
#     F_u_corr_mc = inner(v_mc, u_mc)*dx - inner(v_mc, ut_mc)*dx + dt*g*theta*inner(v_mc, grad(eta_diff_mc))*dx
#     a_u_corr_mc = lhs(F_u_corr_mc)
#     L_u_corr_mc = rhs(F_u_corr_mc)
#
#
#
#     # Assemble matrices
#     A_u_corr_mc = assemble(a_u_corr_mc)
#     for bc_mc in u_bc_object_list_mc: bc_mc.apply(A_u_corr_mc)
#     a_u_corr_solver_mc = LUSolver(A_u_corr_mc)
#     a_u_corr_solver_mc.parameters["reuse_factorization"] = True
#
#     if linear_divergence:
#         A_p_corr_mc = assemble(a_p_corr_mc)
#         for bc_mc in eta_bc_object_list_mc: bc_mc.apply(A_p_corr_mc)
#         a_p_corr_solver_mc = LUSolver(A_p_corr_mc)
#         a_p_corr_solver_mc.parameters["reuse_factorization"] = True
#
#     #######################################################################
#     # Solving the equation starts ...
#     #######################################################################
#
#     # Time stepping
#     timestepcount_mc = 0
#     result_u1 = []
#     result_v1 = []
#     result_eta1 = []
#     result_u1.append(u1_mc(testnodes[0][0], testnodes[0][1])[0])
#     result_v1.append(u1_mc(testnodes[0][0], testnodes[0][1])[1])
#     result_eta1.append(eta1_mc(testnodes[0][0], testnodes[0][1]))
#
#     while float(t_mc - finish_time) < -1e-3:
#
#         # update current time t.==> successfully updated
#         timestepcount_mc += 1
#         t_mc = Constant(t_mc+dt)
#
#         # update force term with time t. ==> successfully updated. ==> weak force implemented here.
#         t_theta_mc = Constant(t_mc-(1-theta)*dt)
#         f_u_expression_object_list_mc.t = t_theta_mc
#         f_v_expression_object_list_mc.t = t_theta_mc
#
#         # update boundary term with time t. ==> successfully updated. ==> strong dirichlet BC implemented here.
#         if uTimeDependent_mc:
#             u_list_expression_mc.t = t_mc
#         if etaTimeDependent_mc:
#             eta_list_expression_mc.t = t_mc
#
#         # update eddy_viscosity
#         if include_les:
#             les_mc.u = u0_mc
#             les_mc.solve()
#
#         # Compute tentative velocity step.
#         A_u_tent_mc = assemble(a_u_tent_mc)
#         b_mc = assemble(L_u_tent_mc)
#         for bc_mc in u_bc_object_list_mc: bc_mc.apply(A_u_tent_mc, b_mc)
#         solve(A_u_tent_mc, ut_mc.vector(), b_mc)
#
#         # Compute pressure correction step.
#         b_mc = assemble(L_p_corr_mc)
#         for bc_mc in eta_bc_object_list_mc: bc_mc.apply(b_mc)
#         if linear_divergence:
#             a_p_corr_solver_mc.solve(eta1_mc.vector(), b_mc)
#         else:
#             A_p_corr_mc = assemble(a_p_corr_mc)
#             for bc_mc in eta_bc_object_list_mc: bc_mc.apply(A_p_corr_mc)
#             solve(A_p_corr_mc, eta1_mc.vector(), b_mc)
#
#         # Compute Velocity correction step.
#         b_mc = assemble(L_u_corr_mc)
#         for bc_mc in u_bc_object_list_mc: bc_mc.apply(b_mc)
#         a_u_corr_solver_mc.solve(u1_mc.vector(), b_mc)
#
#         # Rotate functions for next time_step
#         u00_mc.assign(u0_mc)
#         u0_mc.assign(u1_mc)
#         eta0_mc.assign(eta1_mc)
#
#         # save to a result list
#         result_u1.append(u1_mc(testnodes[0][0], testnodes[0][1])[0])
#         result_v1.append(u1_mc(testnodes[0][0], testnodes[0][1])[1])
#         result_eta1.append(eta1_mc(testnodes[0][0], testnodes[0][1]))
#
#
#     return result_u1, result_v1, result_eta1


# binrandom_u1_true = np.zeros([n_sample, n_sample, ntimestep+1])
# binrandom_v1_true = np.zeros([n_sample, n_sample, ntimestep+1])
# binrandom_eta1_true = np.zeros([n_sample, n_sample, ntimestep+1])
#
# sample_initial_eta_key = initial_eta.keys()
# sample_initial_u_key = initial_u.keys()
# sample_boundary_u_key = boundary_u.keys()
# sample_boundary_eta_key = boundary_eta.keys()
# sample_bathy_key = bathymetry.keys()
#
# sample_initial_eta_value = initial_eta.values()
# sample_initial_u_value = initial_u.values()
# sample_boundary_u_value = boundary_u.values()
# sample_boundary_eta_value = boundary_eta.values()
# sample_bathy_value = bathymetry.values()

# for j in range(n_sample):
#     for k in range(n_sample):
#         # recal sto_viscosity sto_force
#         q0 = samplex[j]
#         q1 = sampley[k]
#         x = sp.Symbol('x')
#         y = sp.Symbol('y')
#         sample_viscosity = str(eval(sto_viscosity)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         sample_friction = str(eval(sto_friction)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         sample_force_x = str(eval(sto_force_x)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         sample_force_y = str(eval(sto_force_y)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         # if sample_bathy.has_key('vary'):
#         #     temp = sample_bathy.get('vary')
#         #     temp1 = str(eval(temp)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     sample_bathy['vary'] = temp1
#         # elif sample_bathy.has_key('class') and sample_bathy.values()[0][0] == 'type1':
#         #     temp = sample_bathy.get('class')
#         #     temp1 = str(eval(temp[3])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     temp2 = str(eval(temp[4])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     sample_bathy['class'][3] = temp1
#         #     sample_bathy['class'][4] = temp2
#         # else:
#         #     pass
#         # if sample_initial_eta.has_key('vary'):
#         #     temp = initial_eta.get('vary')
#         #     temp1 = str(eval(temp)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     sample_initial_eta['vary'] = temp1
#         # if sample_initial_u.has_key('vary'):
#         #     temp = initial_u.get('vary')
#         #     temp1 = str(eval(temp[0])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     temp2 = str(eval(temp[1])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #     sample_initial_u['vary'] = (temp1, temp2)
#         # for i, key in enumerate(sample_boundary_u.keys()):
#         #     value = sample_boundary_u[key]
#         #     if isinstance(value, tuple) and isinstance(value[0], str):
#         #         temp1 = str(eval(value[0])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #         temp2 = str(eval(value[1])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #         sample_boundary_u[key][i] = (temp1, temp2)
#         # for i, key in enumerate(sample_boundary_eta.keys()):
#         #     value = sample_boundary_eta[key]
#         #     if isinstance(value, str) and value != 'freeslipyy' and value != 'freeslipxx' and value != 'noslip':
#         #         temp = str(eval(value)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #         sample_boundary_eta[key][i] = temp
#         sample_bathy = dict(zip(sample_bathy_key, sample_bathy_value))
#         sample_initial_eta = dict(zip(sample_initial_eta_key, sample_initial_eta_value))
#         sample_initial_u = dict(zip(sample_initial_u_key, sample_initial_u_value))
#         sample_boundary_eta = dict(zip(sample_boundary_eta_key, sample_boundary_eta_value))
#         sample_boundary_u = dict(zip(sample_boundary_u_key, sample_boundary_u_value))
#
#         # test2_enable:
#         #sample_initial_eta_value = str(eval(initial_eta.get('vary'))).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
#         #sample_initial_eta = {sample_initial_eta_key[0]: sample_initial_eta_value}
#
#         # test3_enable:
#         # sample_bathy_value = bathymetry.get('class')[0:3]
#         # temp3 = str(eval(bathymetry.get('class')[3]))
#         # temp4 = str(eval(bathymetry.get('class')[4]))
#         # sample_bathy_value.append(temp3)
#         # sample_bathy_value.append(temp4)
#         # sample_bathy = {sample_bathy_key[0]: sample_bathy_value}
#
#         # test4_enable:
#         # none here, all above
#
#         # give u1 and eta1 at each time step here.
#         result_u1_list, result_v1_list, result_eta1_list = func_mc(sample_viscosity, sample_friction, sample_force_x, sample_force_y, sample_bathy, sample_initial_u, sample_initial_eta, sample_boundary_u, sample_boundary_eta)
#         for i in range(ntimestep+1):
#             binrandom_u1_true[j, k, i] = result_u1_list[i]
#             binrandom_v1_true[j, k, i] = result_v1_list[i]
#             binrandom_eta1_true[j, k, i] = result_eta1_list[i]
#
#         print "MC finish: " + str(j) + str(k)

# save MC and PC to file here.


# np.save(output_dir+"binrandom_eta1_true_order_"+str(sto_poly_deg), binrandom_eta1_true)
# np.save(output_dir+"binrandom_u1_true_order_"+str(sto_poly_deg), binrandom_u1_true)
# np.save(output_dir+"binrandom_v1_true_order_"+str(sto_poly_deg), binrandom_v1_true)

# for i in range(ntimestep+1):
#
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     surf = ax.plot_surface(sampleX, sampleY, bin_random_eta1[:, :, i], linewidth=0, antialiased=False)
#     #fig.colorbar(surf, shrink=0.5, aspect=10)
#     ax.scatter(sampleX, sampleY, binrandom_eta1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
#     ax.set_xlim(coefficient[0], coefficient[1])
#     ax.set_ylim(coefficient[2], coefficient[3])
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zlim(zmin, zmax)
#     ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
#     ax.set_xlabel(r'$\xi_1$')
#     ax.set_ylabel(r'$\xi_2$')
#     ax.set_zlabel('ele')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_3Dsurface_eta_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()
#
#     temp = bin_random_eta1[:, :, i]-binrandom_eta1_true[:, :, i]
#     plt.figure()
#     plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#     cbar = plt.colorbar(format='%.2e')
#     cbar.ax.set_title(r'$eta_{pc}-eta_{mc}$')
#     plt.xlim(coefficient[0], coefficient[1])
#     plt.ylim(coefficient[2], coefficient[3])
#     plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     plt.xlabel(r'$\xi_1$')
#     plt.ylabel(r'$\xi_2$')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_2Ddifference_eta_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()
#
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     surf = ax.plot_surface(sampleX, sampleY, bin_random_u1[:, :, i], linewidth=0, antialiased=False)
#     #fig.colorbar(surf, shrink=0.5, aspect=10)
#     ax.scatter(sampleX, sampleY, binrandom_u1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
#     ax.set_xlim(coefficient[0], coefficient[1])
#     ax.set_ylim(coefficient[2], coefficient[3])
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zlim(zmin, zmax)
#     ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
#     ax.set_xlabel(r'$\xi_1$')
#     ax.set_ylabel(r'$\xi_2$')
#     ax.set_zlabel('u')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_3Dsurface_u_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()
#
#     temp = bin_random_u1[:, :, i]-binrandom_u1_true[:, :, i]
#     plt.figure()
#     plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#     cbar = plt.colorbar(format='%.2e')
#     cbar.ax.set_title(r'$u_{pc}-u_{mc}$')
#     plt.xlim(coefficient[0], coefficient[1])
#     plt.ylim(coefficient[2], coefficient[3])
#     plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     plt.xlabel(r'$\xi_1$')
#     plt.ylabel(r'$\xi_2$')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_2Ddifference_u_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()
#
#     fig = plt.figure()
#     ax = fig.gca(projection='3d')
#     surf = ax.plot_surface(sampleX, sampleY, bin_random_v1[:, :, i], linewidth=0, antialiased=False)
#     #fig.colorbar(surf, shrink=0.5, aspect=10)
#     ax.scatter(sampleX, sampleY, binrandom_v1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
#     ax.set_xlim(coefficient[0], coefficient[1])
#     ax.set_ylim(coefficient[2], coefficient[3])
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zlim(zmin, zmax)
#     ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
#     #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
#     ax.set_xlabel(r'$\xi_1$')
#     ax.set_ylabel(r'$\xi_2$')
#     ax.set_zlabel('v')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_3Dsurface_v_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()
#
#     temp = bin_random_v1[:, :, i]-binrandom_v1_true[:, :, i]
#     plt.figure()
#     plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#     cbar = plt.colorbar(format='%.2e')
#     cbar.ax.set_title(r'$v_{pc}-v_{mc}$')
#     plt.xlim(coefficient[0], coefficient[1])
#     plt.ylim(coefficient[2], coefficient[3])
#     plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
#     plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
#     plt.xlabel(r'$\xi_1$')
#     plt.ylabel(r'$\xi_2$')
#     plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*time_step) + "s")
#     #plt.show()
#     plt.savefig(output_dir+"PCMC_2Ddifference_v_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
#     plt.close()