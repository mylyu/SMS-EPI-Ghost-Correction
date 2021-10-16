For a quick start, run Step4_SMS_PEC_SENSE_run.m directly.

To obtain intermediate results:
1. run Step1_MBGC_VcSake.m to get clean single-band images.
2. run Step2_MBGC_calcAndSaveCsm.m to calculate coil sensitivity maps.
3. run Step3_MBGC_CalcAndSavePemUsingCsm.m to separately reconstruct the pos/neg echoes of single-band data (for phase error map calculation).

Note that for 1d LPC on single-band data, LmyGhostCorrection.oneDimLinearCorr_entropy can be used.

-------------- Requirements---------------
1. BART (https://mrirecon.github.io/bart/)
2. Matlab 2016b or above (maybe older versions are ok, but I did not test)
3. ncigt_fil_v2.2_20150119 (for pos_neg_add function)
4. ESPIRiT toolbox (attached)

Let me know (lvmengye@gmail.com) if there is any missing function 