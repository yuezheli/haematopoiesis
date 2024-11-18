function scdODEs!(du,u,p,t)

    # endogenous branch 
    LT0 = u[1]
    ST01 = u[2]
    ST02 = u[3]
    MPP01 = u[4]
    MPP02 = u[5]
    CMP01 = u[6]
    CMP02 = u[7]
    BFUE01 = u[8]
    BFUE02 = u[9]
    CFUE01 = u[10]
    CFUE02 = u[11]
    RET0 = u[12]
    RBC0 = u[13]
    alpha_RET0 = u[14]
    beta_RET0 = u[15]
    gamma_RET0 = u[16]
    delta_RET0 = u[17]
    alphabeta_RET0 = u[18]
    alphagamma_RET0 = u[19]
    alphadelta_RET0 = u[20]
    HbS_RET0 = u[21]
    HbF_RET0 = u[22]
    HbA2_RET0 = u[23]
    alpha_RBC0 = u[24]
    beta_RBC0 = u[25]
    gamma_RBC0 = u[26]
    delta_RBC0 = u[27]
    alphabeta_RBC0 = u[28]
    alphagamma_RBC0 = u[29]
    alphadelta_RBC0 = u[30]
    HbS_RBC0 = u[31]
    HbF_RBC0 = u[32]
    HbA2_RBC0 = u[33]

    # transduced branch,  i = 1,2,3,4
    LT1 = u[34:37]
    ST1 = u[38:41]
    ST2 = u[42:45]
    MPP1 = u[46:49]
    MPP2 = u[50:53]
    CMP1 = u[54:57]
    CMP2 = u[58:61]
    BFUE1 = u[62:65]
    BFUE2 = u[66:69]
    CFUE1 = u[70:73]
    CFUE2 = u[74:77]
    RET = u[78:81]
    RBC = u[82:85]
    alpha_RET = u[86:89]
    beta_RET = u[90:93]
    betanew_RET = u[94:97]
    gamma_RET = u[98:101]
    delta_RET = u[102:105]
    alphabeta_RET = u[106:109]
    alphabetanew_RET = u[110:113]
    alphagamma_RET = u[114:117]
    alphadelta_RET = u[118:121]
    HbS_RET = u[122:125]
    HbA_RET = u[126:129]
    HbF_RET = u[130:133]
    HbA2_RET = u[134:137]
    alpha_RBC = u[138:141]
    beta_RBC = u[142:145]
    betanew_RBC = u[146:149]
    gamma_RBC = u[150:153]
    delta_RBC = u[154:157]
    alphabeta_RBC = u[158:161]
    alphabetanew_RBC = u[162:165]
    alphagamma_RBC = u[166:169]
    alphadelta_RBC = u[170:173]
    HbS_RBC = u[174:177]
    HbA_RBC = u[178:181]
    HbF_RBC = u[182:185]
    HbA2_RBC = u[186:189]

    # create ODE differentiation for later use
#=
    dxdt_LT0 = similar(LT0)
    dxdt_ST01 = similar(ST01)
    dxdt_ST02 = similar(ST02)
    dxdt_MPP01 = similar(MPP01)
    dxdt_MPP02 = similar(MPP02)
    dxdt_CMP01 = similar(CMP01)
    dxdt_CMP02 = similar(CMP02)
    dxdt_BFUE01 = similar(BFUE01)
    dxdt_BFUE02 = similar(BFUE02)
    dxdt_CFUE01 = similar(CFUE01)
    dxdt_CFUE02 = similar(CFUE02)
    dxdt_RET0 = similar(RET0)
    dxdt_RBC0 = similar(RBC0)
    dxdt_alpha_RET0 = similar(alpha_RET0)
    dxdt_beta_RET0 = similar(beta_RET0)
    dxdt_gamma_RET0 = similar(gamma_RET0)
    dxdt_delta_RET0 = similar(delta_RET0)
    dxdt_alphabeta_RET0 = similar(alphabeta_RET0)
    dxdt_alphagamma_RET0 = similar(alphagamma_RET0)
    dxdt_alphadelta_RET0 = similar(alphadelta_RET0)
    dxdt_HbS_RET0 = similar(HbS_RET0)
    dxdt_HbF_RET0 = similar(HbF_RET0)
    dxdt_HbA2_RET0 = similar(HbA2_RET0)
    dxdt_alpha_RBC0 = similar(alpha_RBC0)
    dxdt_beta_RBC0 = similar(beta_RBC0)
    dxdt_gamma_RBC0 = similar(gamma_RBC0)
    dxdt_delta_RBC0 = similar(delta_RBC0)
    dxdt_alphabeta_RBC0 = similar(alphabeta_RBC0)
    dxdt_alphagamma_RBC0 = similar(alphagamma_RBC0)
    dxdt_alphadelta_RBC0 = similar(alphadelta_RBC0)
    dxdt_HbS_RBC0 = similar(HbS_RBC0)
    dxdt_HbF_RBC0 = similar(HbF_RBC0)
    dxdt_HbA2_RBC0 = similar(HbA2_RBC0)
=#
    dxdt_LT1 = similar(LT1)
    dxdt_ST1 = similar(ST1)
    dxdt_ST2 = similar(ST2)
    dxdt_MPP1 = similar(MPP1)
    dxdt_MPP2 = similar(MPP2)
    dxdt_CMP1 = similar(CMP1)
    dxdt_CMP2 = similar(CMP2)
    dxdt_BFUE1 = similar(BFUE1)
    dxdt_BFUE2 = similar(BFUE2)
    dxdt_CFUE1 = similar(CFUE1)
    dxdt_CFUE2 = similar(CFUE2)
    dxdt_RET = similar(RET)
    dxdt_RBC = similar(RBC)
    dxdt_alpha_RET = similar(alpha_RET)
    dxdt_beta_RET = similar(beta_RET)
    dxdt_betanew_RET = similar(betanew_RET)
    dxdt_gamma_RET = similar(gamma_RET)
    dxdt_delta_RET = similar(delta_RET)
    dxdt_alphabeta_RET = similar(alphabeta_RET)
    dxdt_alphabetanew_RET = similar(alphabetanew_RET)
    dxdt_alphagamma_RET = similar(alphagamma_RET)
    dxdt_alphadelta_RET = similar(alphadelta_RET)
    dxdt_HbS_RET = similar(HbS_RET)
    dxdt_HbA_RET = similar(HbA_RET)
    dxdt_HbF_RET = similar(HbF_RET)
    dxdt_HbA2_RET = similar(HbA2_RET)
    dxdt_alpha_RBC = similar(alpha_RBC)
    dxdt_beta_RBC = similar(beta_RBC)
    dxdt_betanew_RBC = similar(betanew_RBC)
    dxdt_gamma_RBC = similar(gamma_RBC)
    dxdt_delta_RBC = similar(delta_RBC)
    dxdt_alphabeta_RBC = similar(alphabeta_RBC)
    dxdt_alphabetanew_RBC = similar(alphabetanew_RBC)
    dxdt_alphagamma_RBC = similar(alphagamma_RBC)
    dxdt_alphadelta_RBC = similar(alphadelta_RBC)
    dxdt_HbS_RBC = similar(HbS_RBC)
    dxdt_HbA_RBC = similar(HbA_RBC)
    dxdt_HbF_RBC = similar(HbF_RBC)
    dxdt_HbA2_RBC = similar(HbA2_RBC)

    # parameters for the model
    # total mean residence time (days);
    # mean residence time in each subcompartment will be computed later;
    tauLT = 100
    tauST = 20
    tauMPP = 2
    tauCMP = 4
    tauBFUE = 7
    tauCFUE = 7
    tauRET = 3
    # tauRBChealthy = 120  # this parameter will be changed from outside
    tauRBCsickle = 12
    tauRBCtrans = 90

    ssLT = 1275 # total capacity for LT-HSC

    # amplification parameters (A.U.);
    aST = 1000
    aMPP = 1000
    aCMP = 16
    aBFUE = 32
    # aCFUE = 32 # see appendix for more information. This implementation has issue

    # replication rate for LT-HSC
    rLT = (1/(2.5*7)); # this is an assumed rate. The rate is missing from the original paper; this number should be greater than 0.01


    bloodvolume = 5  # unit in L

    # hemoglobin related parameters
    VRET = 0.09e-12 # reticulocyte volume; L 
    VRBC = 0.09e-12 # erythrocyte volume; L
    MWHh = 64500 # hemoglobin molecular weight, g/mol or unit in ng/nmol

    # change MWHh, following suggestions in Pittman, 2011
    # https:#www.ncbi.nlm.nih.gov/books/NBK54103/
    # MWHh = 64400

    ksynalpha = 6e-7 # alpha globin synthesis rate; nmol/day/cell
    ratiosynbeta = 0.5 # beta globin/ alpha globin synthesis rate ratio; unitless; ratio assumed based on gene copy number
    ratiosynbetanew = 0.5
    ratiosyngamma = 0.03 # gamma globin/ alpha globin synthesis rate ratio
    ratiosyndelta = 0.04 # delta globin/ alpha globin synthesis rate ratio
    thalfmonomer = 0.25 # free monomer half life; days
    kon = 1e-5*(60*60*24) # bimolecular binding on rate constant; unit nmol-1.day-1 (value adjusted for unit)
    Kdalphabeta_trans = 1e-3 # αβ dimer dissociation constant, nM; this is for either healthy or transgene
    Kdalphabeta_sickle = 1e-2 # αβˢ dimer dissociation constant, nM; this is for sickle gene
    Kdalphagamma = 1e-5 # αγ dimer dissociation constant, nM
    Kdalphadelta = 1e-2  # αδ dimer dissociation constant, nM
    KdHbS = 100 # HbS (α₂β₂ˢ) tetramer dissociation constant, nM
    KdHbA = 100 # HbA (α₂β₂) tetramer dissociation constant, nM
    KdHbF = 100 # HbF (α₂γ₂) tetramer dissociation constant, nM
    KdHbA2 = 100 # HbA2 (α₂δ₂) tetramer dissociation rate, nM


    # these parameters are added for further optimization
    HbA_saturation = 0.74
    HbS_saturation = 0.68
    HbF_saturation = 0.88
    HbA2_saturation = 0


    # LT-HSC dynamics
    kLT2ST = 1/tauLT;
    deltaLT = (rLT - kLT2ST)/ssLT; 

    # ST-HSC dynamics
    tau1ST = (log2(aST) -1)/log2(aST) * tauST;
    tau2ST = tauST - tau1ST;
    krep1ST = (aST/2-1)/ tau1ST;
    krep2ST = 1/ tau2ST;
    k12ST = aST/(2*tau1ST);
    kST2MPP = 2/tau2ST;

    # MPP dynamics
    tau1MPP = (log2(aMPP)-1)/log2(aMPP) * tauMPP;
    tau2MPP = tauMPP - tau1MPP;
    krep1MPP = (aMPP/2-1)/ tau1MPP;
    krep2MPP = 1/ tau2MPP;
    k12MPP = aMPP/(2*tau1MPP);
    kMPP2CMP = 2/tau2MPP;

    # CMP dynamics
    tau1CMP = (log2(aCMP)-1)/log2(aCMP) * tauCMP;
    tau2CMP = tauCMP - tau1CMP;
    krep1CMP = (aCMP/2-1)/ tau1CMP;
    krep2CMP = 1/ tau2CMP;
    k12CMP = aCMP/(2*tau1CMP);
    kCMP2BFUE = 2/tau2CMP;

    # BFU-E dynamics
    tau1BFUE = (log2(aBFUE)-1)/log2(aBFUE) * tauBFUE;
    tau2BFUE = tauBFUE - tau1BFUE;
    krep1BFUE = (aBFUE/2-1)/ tau1BFUE;
    krep2BFUE = 1/ tau2BFUE;
    k12BFUE = aBFUE/(2*tau1BFUE);
    kBFUE2CFUE = 2/tau2BFUE;

    # RET and RBC dynamics
    kRET2RBC = 1/tauRET;
    tauRBCendo = tauRBCsickle;
    kdeathRBCendo = 1/tauRBCendo;
    kdeathRBCtrans = 1/tauRBCtrans;


    # monomer synthesis rate
    ksynbeta = ratiosynbeta * ksynalpha;
    ksyngamma = ratiosyngamma * ksynalpha;
    ksyndelta = ratiosyndelta * ksynalpha;

    # monomer degredation rate
    kdeg = log(2)/ thalfmonomer; # all the monomer half life time is the same, so only 1 parameter is used

    # dimer dissociation rate
    koffalphabetatrans = Kdalphabeta_trans * kon;
    koffalphagamma = Kdalphagamma * kon; 
    koffalphadelta = Kdalphadelta * kon; 
    koffalphabetasickle = Kdalphabeta_sickle * kon; 

    # tetramer dissociation rate
    koffHbA = KdHbA * kon; 
    koffHbS = KdHbS * kon; 
    koffHbF = KdHbF * kon; 
    koffHbA2 = KdHbA2 * kon; 

    # total RET & RBC volume
    VRET0 = VRET * RET0;
    VRETtrans = VRET * RET;
    VRBC0 = VRBC * RBC0;
    VRBCtrans = VRBC * RBC;

    # calculate total hemoglobin concentration
    HbA = HbA_RBC[1] + HbA_RBC[2] + HbA_RBC[3] + HbA_RBC[4]
    HbS = HbS_RBC0 + HbS_RBC[1] + HbS_RBC[2] + HbS_RBC[3] + HbS_RBC[4]
    HbF = HbF_RBC0 + HbF_RBC[1] + HbF_RBC[2] + HbF_RBC[3] + HbF_RBC[4]
    HbA2 = HbA2_RBC0 + HbA2_RBC[1] + HbA2_RBC[2] + HbA2_RBC[3] + HbA2_RBC[4]

    # calculate blood oxygen level
    vO2 = ( HbA_saturation * HbA + HbS_saturation * HbS + HbF_saturation * HbF + HbA2_saturation * HbA2 ) * 1.34;

    aCFUE = 550 * exp( -0.23 * vO2 );

    # CFU-E dynamics
    tau1CFUE = (log2(aCFUE)-1)/log2(aCFUE) * tauCFUE;
    tau2CFUE = tauCFUE - tau1CFUE;
    krep1CFUE = (aCFUE/2-1)/ tau1CFUE;
    krep2CFUE = 1/ tau2CFUE;
    k12CFUE = aCFUE/(2*tau1CFUE);
    kCFUE2RET = 2/tau2CFUE;

    
    # dynamics of endogenous blood cells
    dxdt_LT0 = @. rLT * LT0 - deltaLT * LT0 * (LT0 + LT1[1]) - kLT2ST * LT0; 

    dxdt_ST01 = @. kLT2ST * LT0 + (krep1ST - k12ST) * ST01; 
    dxdt_ST02 = @. k12ST * ST01 + (krep2ST - kST2MPP) * ST02; 

    dxdt_MPP01 = @. kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; 
    dxdt_MPP02 = @. k12MPP * MPP01 + (krep2MPP - kMPP2CMP) * MPP02; 

    dxdt_CMP01 = @. kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01;  
    dxdt_CMP02 = @. k12CMP * CMP01 + (krep2CMP - kCMP2BFUE) * CMP02; 

    dxdt_BFUE01 = @. kCMP2BFUE * CMP02 + (krep1BFUE - k12BFUE) * BFUE01;   
    dxdt_BFUE02 = @. k12BFUE * BFUE01 + (krep2BFUE - kBFUE2CFUE) * BFUE02; 

    dxdt_CFUE01 = @. kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; 
    dxdt_CFUE02 = @. k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;  

    dxdt_RET0 = @. kCFUE2RET * CFUE02 - kRET2RBC * RET0;  

    dxdt_RBC0 = @. kRET2RBC * RET0 - kdeathRBCendo * RBC0;  

    # monomer in endogenous reticulocyte; 
    dxdt_alpha_RET0 = @. ksynalpha * RET0 - (kdeg + kRET2RBC)*alpha_RET0 - (kon*alpha_RET0*beta_RET0 + kon*alpha_RET0*gamma_RET0 + kon*alpha_RET0*delta_RET0) + (koffalphabetasickle*alphabeta_RET0 + koffalphagamma*alphagamma_RET0 + koffalphadelta*alphadelta_RET0)*VRET0;
    dxdt_beta_RET0  = @. ksynbeta  * RET0 - (kdeg + kRET2RBC)*beta_RET0  - kon * alpha_RET0 * beta_RET0  + koffalphabetasickle *VRET0 * alphabeta_RET0  ; 
    dxdt_gamma_RET0 = @. ksyngamma * RET0 - (kdeg + kRET2RBC)*gamma_RET0 - kon * alpha_RET0 * gamma_RET0 + koffalphagamma      *VRET0 * alphagamma_RET0 ;
    dxdt_delta_RET0 = @. ksyndelta * RET0 - (kdeg + kRET2RBC)*delta_RET0 - kon * alpha_RET0 * delta_RET0 + koffalphadelta      *VRET0 * alphadelta_RET0 ; 


    # dimers in endogenous reticulocyte; 
    dxdt_alphabeta_RET0  = @. kon * alpha_RET0 * beta_RET0 -  koffalphabetasickle * alphabeta_RET0  * VRET0 - 2 * kon * alphabeta_RET0  * alphabeta_RET0  + 2 * koffHbS  *VRET0 * HbS_RET0  - kRET2RBC * alphabeta_RET0;
    dxdt_alphagamma_RET0 = @. kon * alpha_RET0 * gamma_RET0 - koffalphagamma      * alphagamma_RET0 * VRET0 - 2 * kon * alphagamma_RET0 * alphagamma_RET0 + 2 * koffHbF  *VRET0 * HbF_RET0  - kRET2RBC * alphagamma_RET0; 
    dxdt_alphadelta_RET0 = @. kon * alpha_RET0 * delta_RET0 - koffalphadelta      * alphadelta_RET0 * VRET0 - 2 * kon * alphadelta_RET0 * alphadelta_RET0 + 2 * koffHbA2 *VRET0 * HbA2_RET0 - kRET2RBC * alphadelta_RET0; 

    # tetramer in endogenous reticulocyte; 
    dxdt_HbS_RET0  = @. kon * alphabeta_RET0  * alphabeta_RET0  - (koffHbS  *VRET0 + kRET2RBC) * HbS_RET0;
    dxdt_HbF_RET0  = @. kon * alphagamma_RET0 * alphagamma_RET0 - (koffHbF  *VRET0 + kRET2RBC) * HbF_RET0;
    dxdt_HbA2_RET0 = @. kon * alphadelta_RET0 * alphadelta_RET0 - (koffHbA2 *VRET0 + kRET2RBC) * HbA2_RET0;

    # monomer in endogenous RBC; 
    dxdt_alpha_RBC0 = @. kRET2RBC * alpha_RET0 - (kdeg + kdeathRBCendo) * alpha_RBC0 - (kon*alpha_RBC0*beta_RBC0 + kon*alpha_RBC0*gamma_RBC0 + kon*alpha_RBC0*delta_RBC0) + (koffalphabetasickle * alphabeta_RBC0 + koffalphagamma * alphagamma_RBC0 + koffalphadelta * alphadelta_RBC0) * VRBC0; 
    dxdt_beta_RBC0 =  @. kRET2RBC * beta_RET0  - (kdeg + kdeathRBCendo) * beta_RBC0  - kon * alpha_RBC0 * beta_RBC0  + koffalphabetasickle *VRBC0 * alphabeta_RBC0  ;
    dxdt_gamma_RBC0 = @. kRET2RBC * gamma_RET0 - (kdeg + kdeathRBCendo) * gamma_RBC0 - kon * alpha_RBC0 * gamma_RBC0 + koffalphagamma      *VRBC0 * alphagamma_RBC0 ;
    dxdt_delta_RBC0 = @. kRET2RBC * delta_RET0 - (kdeg + kdeathRBCendo) * delta_RBC0 - kon * alpha_RBC0 * delta_RBC0 + koffalphadelta      *VRBC0 * alphadelta_RBC0 ;

    # dimers in endogenous RBC
    dxdt_alphabeta_RBC0  = @. kRET2RBC * alphabeta_RET0  - kdeathRBCendo * alphabeta_RBC0  + kon * alpha_RBC0 * beta_RBC0  - koffalphabetasickle *VRBC0 * alphabeta_RBC0  - 2*kon* alphabeta_RBC0 *alphabeta_RBC0  + 2*koffHbS *VRBC0 * HbS_RBC0;
    dxdt_alphagamma_RBC0 = @. kRET2RBC * alphagamma_RET0 - kdeathRBCendo * alphagamma_RBC0 + kon * alpha_RBC0 * gamma_RBC0 - koffalphagamma      *VRBC0 * alphagamma_RBC0 - 2*kon* alphagamma_RBC0*alphagamma_RBC0 + 2*koffHbF *VRBC0 * HbF_RBC0;
    dxdt_alphadelta_RBC0 = @. kRET2RBC * alphadelta_RET0 - kdeathRBCendo * alphadelta_RBC0 + kon * alpha_RBC0 * delta_RBC0 - koffalphadelta      *VRBC0 * alphadelta_RBC0 - 2*kon* alphadelta_RBC0*alphadelta_RBC0 + 2*koffHbA2*VRBC0 * HbA2_RBC0;

    # tetramers in endogenous RBC
    dxdt_HbS_RBC0  = @. kRET2RBC * HbS_RET0  - kdeathRBCendo * HbS_RBC0  + kon*alphabeta_RBC0 *alphabeta_RBC0  - koffHbS * VRBC0 * HbS_RBC0;
    dxdt_HbF_RBC0  = @. kRET2RBC * HbF_RET0  - kdeathRBCendo * HbF_RBC0  + kon*alphagamma_RBC0*alphagamma_RBC0 - koffHbF * VRBC0 * HbF_RBC0; 
    dxdt_HbA2_RBC0 = @. kRET2RBC * HbA2_RET0 - kdeathRBCendo * HbA2_RBC0 + kon*alphadelta_RBC0*alphadelta_RBC0 - koffHbA2* VRBC0 * HbA2_RBC0; 

    # transduced branch cells 
    dxdt_LT1[1:4] = @. rLT * LT1[1:4] - deltaLT * LT1[1:4] * (LT0 + LT1[1:4]) - kLT2ST * LT1[1:4]; 

    dxdt_ST1[1:4] = @. kLT2ST * LT1[1:4] + (krep1ST - k12ST) * ST1[1:4];  
    dxdt_ST2[1:4] = @. k12ST * ST1[1:4] + (krep2ST - kST2MPP) * ST2[1:4]; 

    dxdt_MPP1[1:4] = @. kST2MPP * ST2[1:4] + (krep1MPP - k12MPP) * MPP1[1:4]; 
    dxdt_MPP2[1:4] = @. k12MPP * MPP1[1:4] + (krep2MPP - kMPP2CMP) * MPP2[1:4]; 

    dxdt_CMP1[1:4] = @. kMPP2CMP * MPP2[1:4] + (krep1CMP - k12CMP) * CMP1[1:4];  
    dxdt_CMP2[1:4] = @. k12CMP * CMP1[1:4] + (krep2CMP - kCMP2BFUE) * CMP2[1:4]; 

    dxdt_BFUE1[1:4] = @. kCMP2BFUE * CMP2[1:4] + (krep1BFUE - k12BFUE) * BFUE1[1:4];   
    dxdt_BFUE2[1:4] = @. k12BFUE * BFUE1[1:4] + (krep2BFUE - kBFUE2CFUE) * BFUE2[1:4]; 

    dxdt_CFUE1[1:4] = @. kBFUE2CFUE * BFUE2[1:4] + (krep1CFUE - k12CFUE) * CFUE1[1:4]; 
    dxdt_CFUE2[1:4] = @. k12CFUE * CFUE1[1:4] + (krep2CFUE - kCFUE2RET) * CFUE2[1:4];  

    dxdt_RET[1:4] = @. kCFUE2RET * CFUE2[1:4] - kRET2RBC * RET[1:4]; 

    dxdt_RBC[1:4] = @. kRET2RBC * RET[1:4] - kdeathRBCtrans * RBC[1:4];

    # monomer in transduced RET
    dxdt_alpha_RET[1:4]   = @. ksynalpha                      * RET[1:4] - (kdeg + kRET2RBC)* alpha_RET[1:4]   - (kon*alpha_RET[1:4].*beta_RET[1:4] + kon*alpha_RET[1:4].*betanew_RET[1:4] + kon*alpha_RET[1:4].*gamma_RET[1:4] + kon*alpha_RET[1:4].*delta_RET[1:4]) + (koffalphabetasickle*alphabeta_RET[1:4] + koffalphabetatrans*alphabetanew_RET[1:4] + koffalphagamma*alphagamma_RET[1:4] + koffalphadelta*alphadelta_RET[1:4]).*VRETtrans[1:4];
    dxdt_beta_RET[1:4]    = @. ksynbeta * (1-ratiosynbetanew) * RET[1:4] - (kdeg + kRET2RBC)* beta_RET[1:4]    - kon*alpha_RET[1:4] .* beta_RET[1:4]   + koffalphabetasickle *VRETtrans[1:4] .* alphabeta_RET[1:4]; 
    dxdt_betanew_RET[1:4] = @. ksynbeta * ratiosynbetanew     * RET[1:4] - (kdeg + kRET2RBC)* betanew_RET[1:4] - kon*alpha_RET[1:4] .* betanew_RET[1:4]+ koffalphabetatrans  *VRETtrans[1:4] .* alphabetanew_RET[1:4];
    dxdt_gamma_RET[1:4]   = @. ksyngamma                      * RET[1:4] - (kdeg + kRET2RBC)* gamma_RET[1:4]   - kon*alpha_RET[1:4] .* gamma_RET[1:4]  + koffalphagamma      *VRETtrans[1:4] .* alphagamma_RET[1:4];
    dxdt_delta_RET[1:4]   = @. ksyndelta                      * RET[1:4] - (kdeg + kRET2RBC)* delta_RET[1:4]   - kon*alpha_RET[1:4] .* delta_RET[1:4]  + koffalphadelta      *VRETtrans[1:4] .* alphadelta_RET[1:4];

    # dimers in transduced RET
    dxdt_alphabeta_RET[1:4]    = @. kon * alpha_RET[1:4] .* beta_RET[1:4]    - koffalphabetasickle *VRETtrans[1:4].* alphabeta_RET[1:4]    - 2*kon* alphabeta_RET[1:4]    .^2 + 2*koffHbS *VRETtrans[1:4].* HbS_RET[1:4]  - kRET2RBC * alphabeta_RET[1:4];
    dxdt_alphabetanew_RET[1:4] = @. kon * alpha_RET[1:4] .* betanew_RET[1:4] - koffalphabetatrans  *VRETtrans[1:4].* alphabetanew_RET[1:4] - 2*kon* alphabetanew_RET[1:4] .^2 + 2*koffHbA *VRETtrans[1:4].* HbA_RET[1:4]  - kRET2RBC * alphabetanew_RET[1:4];
    dxdt_alphagamma_RET[1:4]   = @. kon * alpha_RET[1:4] .* gamma_RET[1:4]   - koffalphagamma      *VRETtrans[1:4].* alphagamma_RET[1:4]   - 2*kon* alphagamma_RET[1:4]   .^2 + 2*koffHbF *VRETtrans[1:4].* HbF_RET[1:4]  - kRET2RBC * alphagamma_RET[1:4]; 
    dxdt_alphadelta_RET[1:4]   = @. kon * alpha_RET[1:4] .* delta_RET[1:4]   - koffalphadelta      *VRETtrans[1:4].* alphadelta_RET[1:4]   - 2*kon* alphadelta_RET[1:4]   .^2 + 2*koffHbA2*VRETtrans[1:4].* HbA2_RET[1:4] - kRET2RBC * alphadelta_RET[1:4]; 

    # tetramer in transduced RET
    dxdt_HbS_RET[1:4]  = @. kon * alphabeta_RET[1:4]    .^2 - koffHbS *VRETtrans[1:4] .* HbS_RET[1:4] - kRET2RBC * HbS_RET[1:4];
    dxdt_HbA_RET[1:4]  = @. kon * alphabetanew_RET[1:4] .^2 - koffHbA *VRETtrans[1:4] .* HbA_RET[1:4] - kRET2RBC * HbA_RET[1:4];
    dxdt_HbF_RET[1:4]  = @. kon * alphagamma_RET[1:4]   .^2 - koffHbF *VRETtrans[1:4] .* HbF_RET[1:4] - kRET2RBC * HbF_RET[1:4]; 
    dxdt_HbA2_RET[1:4] = @. kon * alphadelta_RET[1:4]   .^2 - koffHbA2*VRETtrans[1:4] .* HbA2_RET[1:4]- kRET2RBC * HbA2_RET[1:4];


    # monomer in transduced RBC
    dxdt_alpha_RBC[1:4]   = @. kRET2RBC * alpha_RET[1:4]   - (kdeg + kdeathRBCtrans) * alpha_RBC[1:4] - (kon*alpha_RBC[1:4].*beta_RBC[1:4] + kon*alpha_RBC[1:4].*betanew_RBC[1:4] + kon*alpha_RBC[1:4].*gamma_RBC[1:4] + kon*alpha_RBC[1:4].*delta_RBC[1:4]) + (koffalphabetasickle * alphabeta_RBC[1:4] + koffalphabetatrans * alphabetanew_RBC[1:4] + koffalphagamma * alphagamma_RBC[1:4] + koffalphadelta * alphadelta_RBC[1:4]) .*VRBCtrans[1:4]; 
    dxdt_beta_RBC[1:4]    = @. kRET2RBC * beta_RET[1:4]    - (kdeg + kdeathRBCtrans) * beta_RBC[1:4]    - kon * alpha_RBC[1:4] .* beta_RBC[1:4]    + koffalphabetasickle *VRBCtrans[1:4].* alphabeta_RBC[1:4];
    dxdt_betanew_RBC[1:4] = @. kRET2RBC * betanew_RET[1:4] - (kdeg + kdeathRBCtrans) * betanew_RBC[1:4] - kon * alpha_RBC[1:4] .* betanew_RBC[1:4] + koffalphabetatrans  *VRBCtrans[1:4].* alphabetanew_RBC[1:4];
    dxdt_gamma_RBC[1:4]   = @. kRET2RBC * gamma_RET[1:4]   - (kdeg + kdeathRBCtrans) * gamma_RBC[1:4]   - kon * alpha_RBC[1:4] .* gamma_RBC[1:4]   + koffalphagamma      *VRBCtrans[1:4].* alphagamma_RBC[1:4];
    dxdt_delta_RBC[1:4]   = @. kRET2RBC * delta_RET[1:4]   - (kdeg + kdeathRBCtrans) * delta_RBC[1:4]   - kon * alpha_RBC[1:4] .* delta_RBC[1:4]   + koffalphadelta     *VRBCtrans[1:4].* alphadelta_RBC[1:4];


    # dimers in transduced RBC
    dxdt_alphabeta_RBC[1:4]    = @. kRET2RBC * alphabeta_RET[1:4]    - kdeathRBCtrans * alphabeta_RBC[1:4]    + kon * alpha_RBC[1:4] .* beta_RBC[1:4]    - koffalphabetasickle *VRBCtrans[1:4] .* alphabeta_RBC[1:4]    - 2*kon* alphabeta_RBC[1:4]   .^2 + 2*koffHbS *VRBCtrans[1:4] .* HbS_RBC[1:4];
    dxdt_alphabetanew_RBC[1:4] = @. kRET2RBC * alphabetanew_RET[1:4] - kdeathRBCtrans * alphabetanew_RBC[1:4] + kon * alpha_RBC[1:4] .* betanew_RBC[1:4] - koffalphabetatrans  *VRBCtrans[1:4] .* alphabetanew_RBC[1:4] - 2*kon* alphabetanew_RBC[1:4].^2 + 2*koffHbA *VRBCtrans[1:4] .* HbA_RBC[1:4];
    dxdt_alphagamma_RBC[1:4]   = @. kRET2RBC * alphagamma_RET[1:4]   - kdeathRBCtrans * alphagamma_RBC[1:4]   + kon * alpha_RBC[1:4] .* gamma_RBC[1:4]   - koffalphagamma      *VRBCtrans[1:4] .* alphagamma_RBC[1:4]   - 2*kon* alphagamma_RBC[1:4]  .^2 + 2*koffHbF *VRBCtrans[1:4] .* HbF_RBC[1:4];
    dxdt_alphadelta_RBC[1:4]   = @. kRET2RBC * alphadelta_RET[1:4]   - kdeathRBCtrans * alphadelta_RBC[1:4]   + kon * alpha_RBC[1:4] .* delta_RBC[1:4]   - koffalphadelta      *VRBCtrans[1:4] .* alphadelta_RBC[1:4]   - 2*kon* alphadelta_RBC[1:4]  .^2 + 2*koffHbA2*VRBCtrans[1:4] .* HbA2_RBC[1:4];


    # tetramer in transduced RBC
    dxdt_HbS_RBC[1:4]  = @. kRET2RBC * HbS_RET[1:4]  - kdeathRBCtrans * HbS_RBC[1:4]  + kon* alphabeta_RBC[1:4]   .^2 - koffHbS *VRBCtrans[1:4] .* HbS_RBC[1:4];
    dxdt_HbA_RBC[1:4]  = @. kRET2RBC * HbA_RET[1:4]  - kdeathRBCtrans * HbA_RBC[1:4]  + kon* alphabetanew_RBC[1:4].^2 - koffHbA *VRBCtrans[1:4] .* HbA_RBC[1:4];
    dxdt_HbF_RBC[1:4]  = @. kRET2RBC * HbF_RET[1:4]  - kdeathRBCtrans * HbF_RBC[1:4]  + kon* alphagamma_RBC[1:4]  .^2 - koffHbF *VRBCtrans[1:4] .* HbF_RBC[1:4]; 
    dxdt_HbA2_RBC[1:4] = @. kRET2RBC * HbA2_RET[1:4] - kdeathRBCtrans * HbA2_RBC[1:4] + kon* alphadelta_RBC[1:4]  .^2 - koffHbA2*VRBCtrans[1:4] .* HbA2_RBC[1:4]; 

    ## endogenous branch
    du[1] = dxdt_LT0
    du[2] = dxdt_ST01
    du[3] = dxdt_ST02  
    du[4] = dxdt_MPP01  
    du[5] = dxdt_MPP02  
    du[6] = dxdt_CMP01  
    du[7] = dxdt_CMP02  
    du[8] = dxdt_BFUE01  
    du[9] = dxdt_BFUE02  
    du[10] = dxdt_CFUE01  
    du[11] = dxdt_CFUE02  
    du[12] = dxdt_RET0  
    du[13] = dxdt_RBC0  
    du[14] = dxdt_alpha_RET0 
    du[15] = dxdt_beta_RET0  
    du[16] = dxdt_gamma_RET0  
    du[17] = dxdt_delta_RET0  
    du[18] = dxdt_alphabeta_RET0 
    du[19] = dxdt_alphagamma_RET0
    du[20] = dxdt_alphadelta_RET0  
    du[21] = dxdt_HbS_RET0
    du[22] = dxdt_HbF_RET0 
    du[23] = dxdt_HbA2_RET0  
    du[24] = dxdt_alpha_RBC0  
    du[25] = dxdt_beta_RBC0  
    du[26] = dxdt_gamma_RBC0  
    du[27] = dxdt_delta_RBC0  
    du[28] = dxdt_alphabeta_RBC0
    du[29] = dxdt_alphagamma_RBC0  
    du[30] = dxdt_alphadelta_RBC0  
    du[31] = dxdt_HbS_RBC0
    du[32] = dxdt_HbF_RBC0  
    du[33] = dxdt_HbA2_RBC0 

    du[34:37] .= dxdt_LT1 
    du[38:41] .= dxdt_ST1 
    du[42:45] .= dxdt_ST2 
    du[46:49] .= dxdt_MPP1 
    du[50:53] .= dxdt_MPP2 
    du[54:57] .= dxdt_CMP1 
    du[58:61] .= dxdt_CMP2
    du[62:65] .= dxdt_BFUE1 
    du[66:69] .= dxdt_BFUE2 
    du[70:73] .= dxdt_CFUE1 
    du[74:77] .= dxdt_CFUE2 
    du[78:81] .= dxdt_RET 
    du[82:85] .= dxdt_RBC 
    du[86:89] .= dxdt_alpha_RET  
    du[90:93] .= dxdt_beta_RET 
    du[94:97] .= dxdt_betanew_RET  
    du[98:101] .= dxdt_gamma_RET  
    du[102:105] .= dxdt_delta_RET 
    du[106:109] .= dxdt_alphabeta_RET 
    du[110:113] .= dxdt_alphabetanew_RET  
    du[114:117] .= dxdt_alphagamma_RET  
    du[118:121] .= dxdt_alphadelta_RET
    du[122:125] .= dxdt_HbS_RET  
    du[126:129] .= dxdt_HbA_RET  
    du[130:133] .= dxdt_HbF_RET  
    du[134:137] .= dxdt_HbA2_RET  
    du[138:141] .= dxdt_alpha_RBC  
    du[142:145] .= dxdt_beta_RBC 
    du[146:149] .= dxdt_betanew_RBC  
    du[150:153] .= dxdt_gamma_RBC  
    du[154:157] .= dxdt_delta_RBC 
    du[158:161] .= dxdt_alphabeta_RBC 
    du[162:165] .= dxdt_alphabetanew_RBC  
    du[166:169] .= dxdt_alphagamma_RBC  
    du[170:173] .= dxdt_alphadelta_RBC  
    du[174:177] .= dxdt_HbS_RBC  
    du[178:181] .= dxdt_HbA_RBC  
    du[182:185] .= dxdt_HbF_RBC  
    du[186:189] .= dxdt_HbA2_RBC  
end