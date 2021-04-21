#CPA e CPCA
import numpy as np; from scipy import optimize; import numpy.matlib as mtlib

def inicia_const_banco():
    par = {}; #cria o dicionário
    par["eps"] = 1.; par["sig"] = 0.; #CPA original, SRK + associativo, OUTROS VALORES PODEM SER ESCOLHIDOS
    par["es"] = par["eps"] * par["sig"]; par["ems"] = par["eps"] + par["sig"]; par["ems_1"] = par["ems"] - 1.;
    par["inv_s_e"] = 1. / (par["sig"] - par["eps"]); par["es_e_s"] = par["es"] - par["ems"];
    A = par["ems_1"]**3. + 9.*par["ems_1"]**2. + 27.*par["ems"];     B = 27.*(par["ems"]+par["es"]) - 3.*par["ems_1"]**2. - 18.*par["ems_1"];
    C = 3.*(par["ems"]+2.);     raizes = np.roots(np.poly1d([A,B,C,-1.])); raizes = raizes[np.imag(raizes)==0]; raizes = raizes[raizes>0];
    raizes = np.real(raizes)
    if len(raizes) != 1:
        print("raízes =",raizes)
        raise("Mais de uma raiz para omega a partir de eps e sig")
    else:
        par["omg"] = raizes[0];
    par["Zcmod"] = (1. - par["omg"]*par["ems_1"])/3.;
    par["psi"] = 3.*par["Zcmod"]**2. + par["ems"]*par["omg"] - par["es_e_s"]*par["omg"]**2.

    par["R"] = 0.08314462; #R em bar L / K / mol
    par["omg_psi_R"] = par["omg"] / par["psi"] / par["R"];
    par["tolP"] = 1e-5; #tolerância para convergência em pressão (P), sendo o valor máximo para a diferença relativa percentual entre dois valores de P
    par["rhoquasezero"] = 1.e-5; #densidade em mol/L considerada como limite inferior de busca

    #Dados:          b         a0        c1    epsAB    betaAB     Tc       M   esquema
    #em           L/mol   bar(L/mol)^2    -   barL/mol    -        K      g/mol   -
    par["Mdata"] = {
    'n-heptano': [0.125351, 28.49716, 0.93644, 0,      0,        554.17, 100.21, 'nenhum'],
    'água 3B 1': [0.014969, 3.005960, 0.35928, 207.97, 21.3e-3,  647.29, 18.015, '3B'],
    'metanol 3B':[0.0334,   4.5897,   1.0068,  160.70, 34.4e-3,  512.64, 32.042, '3B'],
    'n-pentano': [0.091010, 17.900,   0.81352, 0,      0,        479.44, 72.151, 'nenhum'],
    'etanol 3B': [0.0500,    8.5755,  1.0564,  150.00, 17.3e-3,  513.92, 46.069, '3B'],
    'água 4C':   [0.014,     1.06,    0.55,    173.3,  68.56e-3, 647.29, 18.015, '4C'],
    'etanol 2B': [0.048,     7.49,    0.72,    231.2,  9.42e-3,  513.92, 46.069, '2B'],
    'metanol 2B':[0.0309,    4.0531,  0.4310,  245.91, 16.1e-3,  512.64, 32.042, '2B'],
    'H2S 4C':    [0.029320,  3.1995,  0.82805, 14.64,  0.4862,   373.3,  34.08,  '4C'],
    'H2S 3B':    [0.028950,  3.2991,  0.74067, 37.82,  0.2329,   373.3,  34.08,  '3B'],
    'HAc':       [0.0468,    9.1196,  0.4644,  403.23, 4.5e-3,   591.95, 60.052, '1A']};
    #Referência para os parâmetros
    #Ten Years with the CPA part 1: Ind. Eng. Chem. Res., Vol. 45, No. 14, 2006, 4855-4868
    par["refbanco"] = {
    'n-heptano': 'Ten Years with the CPA part 1 - Table A2',
    'água 3B 1': 'Ten Years with the CPA part 1 - Table 9 and A1',
    'metanol 3B':'Ten Years with the CPA part 1 - Table 6 and A1',
    'n-pentano': 'Ten Years with the CPA part 1 - Table A2',
    'etanol 3B': 'Ten Years with the CPA part 1 - Table 6 and A1',
    'água 4C':   'Braz. J. Chem. Eng., Vol. 35, No. 02, 2018, 363-372 - Table 1',
    'etanol 2B': 'Braz. J. Chem. Eng., Vol. 35, No. 02, 2018, 363-372 - Table 1',
    'metanol 2B':'Ten Years with the CPA part 1 - Table 4',
    'H2S 4C':    'Ind. Eng. Chem. Res. 2006, 45, 7688-7699 - Table 1',
    'H2S 3B':    'Ind. Eng. Chem. Res. 2006, 45, 7688-7699 - Table 1',
    'HAc':       'Ten Years with the CPA part 1 - Table 4 and A1'};

    return par

def iniciaparpuro(par): 
    if par["Subs do Banco"] == True:
        par["b"],par["a0"],par["c1"],par["epsAB"],par["betaAB"],par["Tc"],par["MM"],par["esquema"]= par["Mdata"][par["substância"]];
        par["refpar"] = par["refbanco"][par["substância"]];

    #rhocmp e razão a0/(b * R)
    par["a0_bR"] = par["a0"] / (par["b"] * par["R"]); par["rhocmp"] = 0.999 / par["b"];
    if par["esquema"] != "nenhum":
        #bbetaAB, epsAB_R e mbAB
        par["b_betaAB"] = par["b"] * par["betaAB"]; par["epsAB_R"] = par["epsAB"] / par["R"];
        if par["esquema"] == "1A":
            par["MepsAB"] = [par["epsAB"]]; par["MbetaAB"] = np.array([par["betaAB"]]);
            par["MepsAB_R"] = np.array(par["MepsAB"])/par["R"]; par["Sassoc"] = np.array([1]);
        elif par["esquema"] == "2B":
            par["MepsAB"] = [[0, par["epsAB"]],[par["epsAB"], 0]]; par["MbetaAB"] = np.array([[1e-300, par["betaAB"]],[par["betaAB"], 1e-300]]);
            par["MepsAB_R"] = np.array(par["MepsAB"])/par["R"]; par["Sassoc"] = np.array([1, 1]);
        elif par["esquema"] == "3B":
            par["MepsAB"] = [[0, par["epsAB"]],[par["epsAB"], 0]]; par["MbetaAB"] = np.array([[1e-300, par["betaAB"]],[par["betaAB"], 1e-300]]);
            par["MepsAB_R"] = np.array(par["MepsAB"])/par["R"]; par["Sassoc"] = np.array([2, 1]);
        elif par["esquema"] == "4C":
            par["MepsAB"] = [[0, par["epsAB"]],[par["epsAB"], 0]]; par["MbetaAB"] = np.array([[1e-300, par["betaAB"]],[par["betaAB"], 1e-300]]);
            par["MepsAB_R"] = np.array(par["MepsAB"])/par["R"]; par["Sassoc"] = np.array([2, 2]);

    return par

def converteTcPcwSRK(Tc, Pc, w, par):
    # Original
    # b = 0.08664 * par["R"] * Tc / Pc;
    # a0 = 0.42748 * (par["R"] * Tc)**2 / Pc;
    # c1 = 0.480 + 1.574 * w - 0.176 * w**2;
    # Soave, Chem. Eng. Sci. v. 34, p.225-229 (1978)
    b = 0.08664035 * par["R"] * Tc / Pc;
    a0 = 0.42748025 * (par["R"] * Tc)**2 / Pc;
    c1 = 0.47979 + 1.576 * w - 0.1925 * w**2 + 0.025 * w**3;
    return b, a0, c1

def Xassocpuro(rho,T,par):
    if par["esquema"] == 'nenhum':
        X = 1; eta = 0.; Delta = 0; gref = 0.; rhoDelta = 0;
    else:
        eta = par["b"] * rho / 4.;
        if par["g"] == "CPA":
            gref = (1. - eta/2.) / (1. - eta)**3.;
        elif par["g"] == "sCPA":
            gref = 1. / (1. - 1.9*eta);
        # Delta = gref * (np.exp(par["epsAB_R"]/T)-1.) * par["b_betaAB"];
        Delta = gref * (np.exp(par["epsAB_R"]/T + np.log(par["b_betaAB"])) - par["b_betaAB"]);
        rhoDelta = rho * Delta;
        if rhoDelta < 1.e-15:
            if par["esquema"] == '1A':
                X = np.array([1.]);
            else:
                X = np.array([1., 1.]);
        elif rhoDelta < 1.e300:
            if par["esquema"] == '1A':
                XA = (-1.+np.sqrt(1+4*rhoDelta))/(2*rhoDelta); X = np.array([XA]);
            elif par["esquema"] == '2B':
                XA = (-1.+np.sqrt(1+4*rhoDelta))/(2*rhoDelta); X = np.array([XA, XA]); #XB = XA
            elif par["esquema"] == '3B':
                XA = (-1.+rhoDelta+np.sqrt((1.+rhoDelta)**2.+4.*rhoDelta))/(4.*rhoDelta); X = np.array([XA, 2.*XA-1.]); #XB = XA (1 tipo, ocorrência 2), XC = 2*XA - 1
            elif par["esquema"] == '4C':
                XA = (-1.+np.sqrt(1+8*rhoDelta))/(4*rhoDelta); X = np.array([XA, XA]); #XB = XA (1 tipo, ocorrência 2), XC = XD (outro sítio, ocorrência 2) = XA
        else:
            # raizrhoDelta = np.sqrt(rho * par["b_betaAB"] * gref) * np.exp(par["epsAB_R"]/(2.*T));
            raizrhoDelta = np.sqrt(rho * gref) * np.exp( ( par["epsAB_R"]/T + np.log(par["b_betaAB"]) )/2. );
            if par["esquema"] == '1A':
                XA = 1./raizrhoDelta; X = np.array([XA]);
            elif par["esquema"] == '2B':
                XA =1./raizrhoDelta; X = np.array([XA, XA]); #XB = XA
            elif par["esquema"] == '3B':
                XC = 3./2./raizrhoDelta**2; XA = (XC+1.)/2.; X = np.array([XA, XA, 2.*XA-1.]); #XB = XA, XC = 2*XA - 1
            elif par["esquema"] == '4C':
                XA = 1./np.sqrt(2.)/raizrhoDelta; X = np.array([XA, XA]); #XB = XA (1 tipo, ocorrência 2), XC = XD (outro sítio, ocorrência 2) = XA
    return X, eta, Delta, gref, rhoDelta

def alfaT(T, par):
    Tr = T / par["Tc"]
    return (1 + par["c1"] * (1 - np.sqrt(Tr)))**2

def Zpuro(rho, T, par):
    if par["esquema"] != 'nenhum':
        X, eta, Delta, gref, rhoDelta = Xassocpuro(rho,T,par);
        if par["g"] == "CPA":
            dgrefdeta = (2.5 - eta) / (1. - eta)**4.;
        elif par["g"] == "sCPA":
            dgrefdeta = 1.9 * gref**2.;
        rhodlngrefdrho = eta * dgrefdeta / gref;
        Zassoc = -(1. + rhodlngrefdrho) * par["Sassoc"] @ (1. - X) / 2.;
    else:
        Zassoc = 0.

    brho = par["b"] * rho; q = par["a0_bR"] * alfaT(T, par) / T;
    Zcub = 1. / (1. - brho) - q * brho / (1. + par["eps"] * brho) / (1 + par["sig"] * brho);

    return Zcub + Zassoc

def Bpuro(T, par):
    if par["esquema"] == 'nenhum':
        Bassoc = 0.
    elif par["esquema"] == '1A':
        # Delta_rho0 = (np.exp(par["epsAB_R"]/T) - 1.) * par["betaAB"];
        Delta_rho0 = np.exp( par["epsAB_R"]/T + np.log(par["betaAB"]) ) - par["betaAB"];
        Bassoc = -Delta_rho0 / 2.;
    else:
        # Delta_rho0 = (np.exp(par["MepsAB_R"]/T) - 1.) * par["MbetaAB"];
        Delta_rho0 = np.exp( par["MepsAB_R"]/T + np.log(par["MbetaAB"]) ) - par["MbetaAB"];
        Bassoc = -np.array(par["Sassoc"]) @ Delta_rho0 @ np.array(par["Sassoc"]) / 2.;

    Bcub = par["b"] * (1. - par["a0_bR"] * alfaT(T, par) / T);
    return Bcub + Bassoc

def arespuro(rho, T, par): #Helmholtz/nRT residual "da física" (referência é o gás ideal na mesma rho, T e x e não P, T e x)
    if par["esquema"] != 'nenhum':
        X, eta, Delta, gref, rhoDelta = Xassocpuro(rho,T,par);
        aassoc = np.sum([par["Sassoc"]] @ (np.log(X) + (1. - X) / 2.));
    else:
        aassoc = 0.

    brho = par["b"]*rho; q = par["a0_bR"] * alfaT(T, par) / T;
    I = par["inv_s_e"] * np.log( (1 + par["sig"] * brho) / (1 + par["eps"] * brho) );
    arescub = -np.log(1 - brho) - q * I;

    return arescub + aassoc

def kTalphapuro(rho, T, par):
    # Derivadas da energia de Helmholtz residual
    # rho é a densidade molar, mol/L, T a temperatura, K, x o vetor de frações molares
    # kT a compressibilidade isotérmica, 1/bar, alfa o coeficiente de expansão térmica, 1/K
    # Cp é a capacidade calorífica a pressão constante, J/K/mol
    # u é a velocidade do som, m/s
    delta = 0.00001; #tamanho do passo de derivação
    up = 1. + delta; um = 1. - delta; dd = 2.*delta; dp = 1. + dd; dm = 1. - dd; RT = par["R"] * T;
    Z = Zpuro(rho, T, par); dzdelta = 12.*delta;
    Zp0 = Zpuro(rho*up, T, par); Zm0 = Zpuro(rho*um, T, par); Zdp0 = Zpuro(rho*dp, T, par); Zdm0 = Zpuro(rho*dm, T, par);
    Z0p = Zpuro(rho, T*up, par); Z0m = Zpuro(rho, T*um, par); Z0dp = Zpuro(rho, T*dp, par); Z0dm = Zpuro(rho, T*dm, par);
    kT = 1. / rho / RT / (Z + (-Zdp0 + 8.*Zp0 - 8.*Zm0 + Zdm0) / dzdelta);
    alpha = kT * par["R"] * rho * (Z + (-Z0dp + 8.*Z0p - 8.*Z0m + Z0dm) / dzdelta);

    return kT, alpha, Z

def upuro(rho, T, alpha, kT, Cp, par):
    aux = kT*rho*Cp - T*alpha**2.;
    if aux < 0:
        sai = "u imaginário";
    else:
        sai = np.sqrt( 1.e5 * Cp / ( par["MM"]*aux ) ); #1e5 entra para a conversão de unidades
    return sai

def Topliss(T, fase, par):
    B = Bpuro(T, par);
    if B > 0:
        curva = "C"; rhomin = par["rhoquasezero"]; rhomax = par["rhocmp"]; Pinf = 0.; rhoinf = 0.
    else:
        sol = optimize.root_scalar(fobjraizinf,args=(T, par),method='brentq',bracket=[par["rhoquasezero"],par["rhocmp"]]);
        rhoinf = sol.root; Zinf, rhodZdrhoinf = fobjraizmaxmin_(rhoinf, T, par);
        Pinf = Zinf * rhoinf * par["R"] * T;
        if Zinf + rhodZdrhoinf > 0:
            curva = "B"; rhomin = par["rhoquasezero"]; rhomax = par["rhocmp"]; #pode ser atualiza se P != 0
            if fase == "ambas":
                raise("Curva tipo B, não é possível haver duas fases");
        else:
            curva = "A";
            if fase == "líquido":
                sol = optimize.root_scalar(fobjraizmaxmin,args=(T, par),method='brentq',bracket=[rhoinf,par["rhocmp"]]);
                rhomin = sol.root; rhomax = 0.;
            elif fase == "vapor":
                sol = optimize.root_scalar(fobjraizmaxmin,args=(T, par),method='brentq',bracket=[par["rhoquasezero"],rhoinf]);
                rhomax = sol.root; rhomin = 0;
            elif fase == "ambas":
                sol = optimize.root_scalar(fobjraizmaxmin,args=(T, par),method='brentq',bracket=[rhoinf,par["rhocmp"]]);
                rhomin = sol.root;
                sol = optimize.root_scalar(fobjraizmaxmin,args=(T, par),method='brentq',bracket=[par["rhoquasezero"],rhoinf]);
                rhomax = sol.root;
            else:
                raise("fase deve ser líquido ou vapor ou ambas")
    return curva, rhomin, rhomax, rhoinf, Pinf

def fobjraizinf(rho, T, par):
    if rho <= par["rhoquasezero"]:
        fobj = 2.*Bpuro(T, par);
    else:
        delta = 0.001; #tamanho do passo de derivação
        up = 1. + delta; um = 1. - delta; dd = 2.*delta; dp = 1. + dd; dm = 1. - dd;
        Z0 = Zpuro(rho, T, par); deltarho = delta * rho;
        Zp = Zpuro(rho*up, T, par); Zm = Zpuro(rho*um, T, par); Zdp = Zpuro(rho*dp, T, par); Zdm = Zpuro(rho*dm, T, par);
        dZdrho = (-Zdp + 8.*Zp - 8.*Zm + Zdm) / 12. / deltarho;
        rhod2Zdrho2 = rho * (-Zdp + 16.*Zp - 30.*Z0 + 16.*Zm - Zdm) / 12. / deltarho**2.
        fobj = 2.*dZdrho + rhod2Zdrho2;
    return fobj;

def fobjraizmaxmin(rho, T, par):
    Z, rhodZdrho = fobjraizmaxmin_(rho, T, par);
    return Z + rhodZdrho;

def fobjraizmaxmin_(rho, T, par):
    if rho <= par["rhoquasezero"]:
        Z0 = 1.; rhodZdrho = 0.
    else:
        delta = 0.001; #tamanho do passo de derivação
        up = 1. + delta; um = 1. - delta; dd = 2.*delta; dp = 1. + dd; dm = 1. - dd;
        Z0 = Zpuro(rho, T, par);
        Zp = Zpuro(rho*up, T, par); Zm = Zpuro(rho*um, T, par); Zdp = Zpuro(rho*dp, T, par); Zdm = Zpuro(rho*dm, T, par);
        rhodZdrho = (-Zdp + 8.*Zp - 8.*Zm + Zdm) / 12. / delta
    return Z0, rhodZdrho;

def fobjraizrho(rho, T, P, fase, par):
    return Zpuro(rho, T, par) * rho - P / par["R"] / T;

def raizrho(T, P, fase, par):
    curva, rhomin, rhomax, rhoinf, Pinf = Topliss(T, fase, par);
    if fase == "líquido" or fase == "ambas":
        Pmin = Zpuro(rhomin, T, par)*rhomin*par["R"]*T
        if P < Pmin:
            # print("Não é possível encontrar raiz de líquido, Pmin = ",Pmin);
            sai = "não achou raiz";
            return sai
    if fase == "vapor" or fase == "ambas":
        Pmax = Zpuro(rhomax, T, par)*rhomax*par["R"]*T;
        if P > Pmax:
            print("Não é possível encontrar raiz de vapor, Pmax = ",Pmax); sai = "não achou raiz";
            return sai
    if curva == "B":
        if P > Pinf:
            liminf = rhoinf;
        else:
            limsup = rhoinf;
    if fase == "líquido":
        # rho0 = par["rhocmp"] / 1.3;
        liminf = rhomin; limsup = par["rhocmp"];
        sol = optimize.root_scalar(fobjraizrho,args=(T, P, fase, par),method='brentq',bracket=[liminf,limsup]);
        sai = sol.root;
    elif fase == "vapor":
        # rho0 = P/(par["R"]*T);
        liminf = par["rhoquasezero"]; limsup = rhomax;
        sol = optimize.root_scalar(fobjraizrho,args=(T, P, fase, par),method='brentq',bracket=[liminf,limsup]);
        sai = sol.root;
    elif fase == "ambas":
        # rho0 = par["rhocmp"] / 1.3;
        liminf = rhomin; limsup = par["rhocmp"];
        solL = optimize.root_scalar(fobjraizrho,args=(T, P, fase, par),method='brentq',bracket=[liminf,limsup]);
        liminf = par["rhoquasezero"]; limsup = rhomax;
        solV = optimize.root_scalar(fobjraizrho,args=(T, P, fase, par),method='brentq',bracket=[liminf,limsup]);
        sai = [solL.root, solV.root];
    else:
        sai = "Fase inválida em raizrho"; #print(sai);
        return sai
    return sai

def rhoupuro(T, P, Cp, fase, par):
    rho = raizrho(T, P, fase, par);
    if type(rho) == str:
        v = "não calculado";
    else:
        kT, alpha, Z = kTalphapuro(rho, T, par); v = upuro(rho, T, alpha, kT, Cp, par);
    return rho, v

def fobjPalphakTpuro_(v, tipo, escala, par, Dadosexp):
    Texp, rhoexp, Pexp, alphaexp, kTexp = Dadosexp; par["Subs do Banco"] = False;
    if np.ndim(v) == 1:
        nparam = 1; v = [v];
    else:
        nparam = np.shape(v)[0];
    vfobj = [];
    for i in range(nparam):
        vi = v[i];
        if tipo == "3p": #3 parâmetros estimados: Tc, Pc e omega ou a0, b e c1
            par["b"], par["a0"], par["c1"] = abs(np.array(vi)) * escala; par["esquema"]="nenhum";
            par["c1"] = vi[2] * escala[2];
            par["Tc"] = par["omg_psi_R"]*par["a0"]/par["b"];
        elif tipo == "5p": #5 parâmetros estimados: a0, b, c1, epsAB, betaAB
            par["b"], par["a0"], par["c1"], par["epsAB"] = abs(np.array(vi[0:-1])) * escala[0:-1];
            par["c1"] = vi[2] * escala[2];
            auxbetaAB = vi[-1] * escala[-1]; par["betaAB"] = par["betaABmax"] / (1. + np.exp(auxbetaAB));
            par["Tc"] = par["omg_psi_R"]*par["a0"]/par["b"];
        par = iniciaparpuro(par);

        Pcalc = []; alphacalc = []; kTcalc = [];
        for j in range(len(Texp)):
            kT, alpha, Z = kTalphapuro(rhoexp[j], Texp[j], par);
            Pcalc.append(Z*rhoexp[j]*par["R"]*Texp[j]); alphacalc.append(alpha*1e3); kTcalc.append(kT*1.e4);
        fobj = np.sum(((Pcalc-Pexp)/Pexp)**2. + ((kTcalc-kTexp)/kTexp)**2. )#+ ((alphacalc-alphaexp)/alphaexp)**2.)
        vfobj.append(fobj);
    return vfobj, Pcalc, alphacalc, kTcalc

def fobjPalphakTpuro(v, tipo, escala, par, Dadosexp):
    vfobj, Pcalc, alphacalc, kTcalc = fobjPalphakTpuro_(v, tipo, escala, par, Dadosexp)
    return np.squeeze(vfobj)

def fobjrhoupuro_(v, tipo, escala, par, Dadosexp):
    Texp, Pexp, rhoexp, uexp, Cpexp, dprhoexp, dpuexp = Dadosexp; par["Subs do Banco"] = False;
    if np.ndim(v) == 1:
        nparam = 1; v = [v];
    else:
        nparam = np.shape(v)[0];
    vfobj = [];
    for i in range(nparam):
        vi = v[i];
        if tipo == "3p": #3 parâmetros estimados: Tc, Pc e omega ou a0, b e c1
            par["b"], par["a0"], par["c1"] = abs(np.array(vi)) * escala; par["esquema"]="nenhum";
            par["c1"] = vi[2] * escala[2];
            par["Tc"] = par["omg_psi_R"]*par["a0"]/par["b"];
        elif tipo == "5p": #5 parâmetros estimados: a0, b, c1, epsAB, betaAB
            par["b"], par["a0"], par["c1"], par["epsAB"] = abs(np.array(vi[0:-1])) * escala[0:-1];
            par["c1"] = vi[2] * escala[2];
            auxbetaAB = vi[-1] * escala[-1]; par["betaAB"] = par["betaABmax"] / (1. + np.exp(auxbetaAB));
            par["Tc"] = par["omg_psi_R"]*par["a0"]/par["b"];
        par = iniciaparpuro(par);

        rhocalc = []; ucalc = [];
        for j in range(len(Texp)):
            rho, u = rhoupuro(Texp[j], Pexp[j], Cpexp[j], "líquido", par);
            if type(rho) == str or type(u) == str:
                rhocalc.append(1e10); ucalc.append(1e10);
            else:
                rhocalc.append(rho); ucalc.append(u);
        fobj = np.sum(((rhocalc-rhoexp)/dprhoexp)**2. + ((ucalc-uexp)/dpuexp)**2. );
        vfobj.append(fobj);
    return vfobj, rhocalc, ucalc

def fobjrhoupuro(v, tipo, escala, par, Dadosexp):
    vfobj, rhocalc, ucalc = fobjrhoupuro_(v, tipo, escala, par, Dadosexp)
    return np.squeeze(vfobj)

