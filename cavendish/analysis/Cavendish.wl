(* ::Package:: *)

(* ::Text:: *)
(*This Notebook contains portions of analysis for Advanced Experimental physics Cavendish Experiment.*)
(*The first part calculates the mean and standard deviation of lengths we measured in the lab.*)


(* ::Input:: *)
(*Smeasurements={0.869, 0.869, 0.87, 0.87, 0.87, 0.867};*)
(*Smean=Mean[Smeasurements];*)
(*SstdDev=StandardDeviation[Smeasurements];*)
(*S = Around[Smean, SstdDev]*)
(*Lmeasurements= {6.653,6.624};*)
(*Lmean = Mean[Lmeasurements];*)
(*LstdDev = StandardDeviation[Lmeasurements];*)
(*L = Around[Lmean,LstdDev]*)
(*b1measurements = {1.382, 1.383, 1.378, 1.373, 1.373, 1.374, 1.379, 1.377, 1.375};*)
(*b1measurementsmeters = 0.0254*b1measurements;*)
(*b1mean = Mean[b1measurementsmeters];*)
(*b1stdDev = StandardDeviation[b1measurementsmeters];*)
(*b2measurements = {1.392,1.384,1.392,1.383,1.378,1.377,1.386,1.38,1.381};*)
(*b2measurementsmeters= 0.0254*b2measurements;*)
(*b2mean = Mean[b2measurementsmeters];*)
(*b2stdDev = StandardDeviation[b2measurementsmeters];*)
(*b1 = Around[b1mean, b1stdDev];*)
(*b2 = Around[b2mean, b2stdDev];*)
(*Qmeasurements = {1.382, 1.383, 1.378, 1.373, 1.373, 1.374, 1.379, 1.377, 1.375, 1.392,1.384,1.392,1.383,1.378,1.377,1.386,1.38,1.381};*)
(*Qmeasurementsmeters = 0.0254*Qmeasurements;*)
(*Qmean = Mean[Qmeasurementsmeters];*)
(*QstdDev = StandardDeviation[Qmeasurementsmeters];*)
(*Q = Around[Qmean, QstdDev];*)
(*b = Q/2 + 0.0554/2*)


(* ::Text:: *)
(**)
(*Now that we have found these values, we use them in our python code to run fits and extract parameters from the fit. Once we have those parameters, we can plug them into our notebook along with necessary values given to us from the pendulum manufacturer. *)


(* ::Input:: *)
(*mb = Around[0.00717, 0];*)
(*lb = Around[0.149,0];*)
(*wb = Around[0.0127, 0];*)
(*d=Around[0.0667, 0];*)
(*M = Around[1.508, 0];*)
(*m = Around[0.0145, 0];*)
(*dss = Around[0.0135, 0];*)
(*b = Around[0.045231, 0.0000];*)
(*int = 0.14144;*)


(* ::Input:: *)
(*thetaA = Around[0.00935, 0.00003];*)
(*thetaB = Around[0.00739,0.00003];*)
(*thetaC= Around[0.00946, 0.00003];*)
(*thetaD = Around[0.00739,0.00003];*)
(*tA = Around[241.51749, 0.37140];*)
(*tB = Around[240.41316,0.90445];*)
(*tC = Around[241.35471,0.33248];*)
(*tD = Around[241.42847,0.75519];*)
(*t1 = Mean[{tA,tB}]*)
(*t2 = Mean[{tC, tD}]*)
(**)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Ibeam[m_, l_, w_] :=m * (1/12)(l^2+w^2)*)
(*ISpheres[m_, dss_, d_] := 2*m((2/5)*(dss/2)^2 + d^2)*)
(*ITotal[Ib_,Is_] := Ib + Is*)
(*Lambda[I_, T_]:= (4*(Pi^2)*I)/(T^2)*)
(*DeltaTheta[thetaA_, thetaB_]:= (thetaA-thetaB)*)
(*DeltaTau[lambda_, deltheta_]:= lambda*deltheta*)
(*GFirst[deltau_, b_, d_, M_, m_]:= (deltau*(b^2))/(4*m*M*d)*)
(*Theta1[b_, d_]:= ArcTan[b/(2d)]*)
(*R1[b_,d_]:=Sqrt[b^2+(2d)^2]*)
(*GSecond[deltau_,b_, d_, M_, m_, theta1_, r_]:= deltau/((4*M*m*d)*((1/b^2)-(Sin[theta1]/r^2)))*)
(*GThird[deltau_,b_,d_,M_,m_,theta1_,r_,mb_,lb_]:=deltau/((4*M*m*d)*((1/b^2)-(Sin[theta1]/r^2))+int)*)
(**)


(* ::Input:: *)
(*ibeam = Ibeam[mb,lb,wb];*)
(*ispheres = ISpheres[m,dss,d];*)
(*itotal = ITotal[ibeam,ispheres];*)
(*lambda = Lambda[itotal, t1];*)
(*theta1 = Theta1[b,d];*)
(*r1 = R1[b,d];*)
(*deltheta = DeltaTheta[thetaA, thetaB];*)
(*deltau = DeltaTau[lambda,deltheta];*)
(*gfirst1 =GFirst[deltau,b,d,M,m] *)
(**)


(* ::Text:: *)
(*Second order calculation:*)


(* ::Input:: *)
(*gsecond1=GSecond[deltau,b,d,M,m,theta1,r1]*)


(* ::Input:: *)
(*gthird1 = GThird[deltau,b,d,M,m,theta1,r1, mb, lb]*)


(* ::Input:: *)
(*lambda = Lambda[itotal, t2];*)
(*deltheta = DeltaTheta[thetaC, thetaD];*)
(*deltau = DeltaTau[lambda,deltheta];*)
(*gfirst2 =GFirst[deltau,b,d,M,m] *)


(* ::Input:: *)
(*gsecond2=GSecond[deltau,b,d,M,m,theta1,r1]*)


(* ::Input:: *)
(*gthird2 = GThird[deltau,b,d,M,m,theta1,r1, mb, lb]*)


(* ::Text:: *)
(*Now that we have values of G for each approximation for two separate trials, we can average the trials to get our measurement of G:*)


(* ::Input:: *)
(*G1st = Mean[{gfirst1, gfirst2}]*)


(* ::Input:: *)
(*G2nd = Mean [{gsecond1, gsecond2}]*)


(* ::Input:: *)
(*G3rd = Mean[{gthird1, gthird2}]*)
