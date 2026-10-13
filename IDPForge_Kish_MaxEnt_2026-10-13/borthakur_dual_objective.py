def func_sigma_reg_sigma_md(l, traj_r, obs_exp_r, weight_bias, sigma_reg, sigma_md):
        l = np.array(l)  # ensure array
        weight_all = +weight_bias  # copy
        weight_all -= np.dot(l, traj_r)  # maxent correction
        shift_all = np.max(weight_all)  # shift to avoid overflow
        weight_all = np.exp(weight_all - shift_all)  # compute weights
        weight_0 = +weight_bias  # copy
        shift_0 = np.max(weight_0)  # shift to avoid overflow
        weight_0 = np.exp(weight_0 - shift_0)  # compute weights

        # Gamma function in maxent:
        # Shifts to avoid overflows
        f = np.log(np.sum(weight_all) / np.sum(weight_0)) + shift_all - shift_0 + np.dot(l, obs_exp_r)

        # derivative of Gamma function:
        der = obs_exp_r - np.dot(traj_r, weight_all) / np.sum(weight_all)  # derivative with respect to l

        f += 0.5 * np.sum((sigma_reg * sigma_reg) * l ** 2 + sigma_md * sigma_md * l ** 2)
        der += sigma_reg * sigma_reg * l + sigma_md * sigma_md * l
        return (f, der)
