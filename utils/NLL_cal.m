function NLL = NLL_cal(I_est,I_obs)

NLL = sum(sum(I_est)) - sum(sum(I_obs.*log(I_est)));

end