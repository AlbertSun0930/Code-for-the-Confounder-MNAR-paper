# Code-for-the-Confounder-MNAR-paper
This is the computer programs and data for 
"Identification and Estimation of Causal Effects with Confounders Missing Not at Random” 
by Jian Sun and Bo Fu.

################

In the 'Simulation' folder:

'1_para_ContinuousY.R' contains the codes in Section 11.1 of the supplementary material for illustrating the WEE method in the continuous outcome scenario when Assumptions (2.1)-(2.5) are satisfied. 

'2_para_BinaryY.R' contains the codes in Section 11.1 of the supplementary material for illustrating the WEE method in the binary outcome scenario when Assumptions (2.1)-(2.5) are satisfied. 

'3_para_Sen.R' contains the codes in Section 11.2 of the supplementary material for investigating the robustness of the WEE method when the treatment-independent missing assumption is violated.

'4_ATE_est.R' contains the codes in Section 4.2 of the main text and Section 11.3 of the supplementary material for illustrating the WEE-based average treatment effect estimators when Assumptions (2.1)-(2.5) are satisfied. 

'5_ATE_Sen.R' contains the codes in Section 11.4 of the supplementary material for investigating the robustness of the WEE-based causal effect estimators when the treatment-independent missing assumption is violated.

################

In the 'Real data analysis' folder:

'RDA_Code.R' contains codes in Section 5 of the main text for illustrating the average treatment effect of marital status on mental health with NHANES data. 
