import sys
import math
from pathlib import Path

#def bppt_bpscore(cbps):
#    cbpsc = 1
#    cbpsN = list(cbps)
#    for i in range(0,len(cbpsN)):
#        if i in PWMBP:
#            cbpsc = cbpsc*PWMBP[i][cbpsN[i]]
#        else:
#            print("Error: the input bps is longer than the PWM")
#    return(cbpsc)
#
#def bppt_get_bpscore(l,basebp):
#    nn = ['A','C','G','T'] # nucleotides
#    NNS = ['A','C','G','T'] #  nucleotides combination with fixed length
#    NN = [] #  nucleotides combination with fixed length
#    
#    for i in range(2,l+1):
#        for nns in NNS:
#            for n in nn:
#                newN = nns+n
#                NN.append(newN)
#        NNS = NN 
#        NN = []
#
#    basebpsc = bppt_bpscore(basebp)
#    cBPSC = {}
#    for nn in NNS:
#        cBPSC[nn] = bppt_bpscore(nn)/basebpsc
#    
#    return(cBPSC)
#
#def bppt_get_AGEZ(seq,offset=12,maxL=-1): # get the AGEZ 
#    sL = len(seq)
#    ss = seq[0:-offset].split("AG")
#    if maxL > 0:
#        pAG = sL - maxL
#    else:
#        pAG = sL - (offset+len(ss[-1])+14)
#    if pAG < 0:
#        pAG = 0
#    return(pAG)
#
#def bppt_get_pptsc(pptS,lppt,l_max_ppt):
#    end = len(pptS)-lppt
#    if end > l_max_ppt: # if the ppt sequence is larger than l_max_ppt, only check the l_max_ppt
#        end = l_max_ppt
#    pptbasesc = (l_max_ppt-lppt+1)*PPTS["TTTTTTTT"]
#    pptSC = 0
#    for i in range(0,end):
#        cppts = pptS[i:i+lppt]
#        pptSC += PPTS[cppts]
#    if pptbasesc == 0:
#        return(pptSC)
#    else:
#        return((pptSC/pptbasesc))
#
#def bppt_dis_pro(pAG,offset=22,interval=4):#4 best
#    if abs(pAG-offset) > 700:
#        return(0)
#    else:
#        return(1/(math.exp((abs(pAG-offset))/interval)+1))#best
#
#def bppt_get_BPPTsc(seq,maxL,baseppt): # get the candidate bps and ppt and their score
#    
#    lmotif = 7
#    lppt = 8
#    l_max_ppt = 20
#    pstart = bppt_get_AGEZ(seq=seq,offset=12,maxL=maxL)
#    sL = len(seq)
#
#    totbpsc = 0
#    totpptsc = 0
#    totsc = 0
#    npos = 0
#
#    orinp = []
#    orsc = []
#
#    for ipos in range(pstart,sL-3-lmotif): # the BPS could close to the 3' end with only one nucleotide distance 
#    #for ipos in range(pstart,sL-14-lmotif): # the BPS could close to the 3' end with only one nucleotide distance 
#        pAG = sL - ipos - 5
#        bpS = seq[ipos:ipos+lmotif] # bps sequence
#        bpSC = cBPSC[bpS] # bps sequence score
#
#        pptSC = 0
#        dis3 = pAG - 1 # the distance of BPS last nucleotide to the 3' end
#        if dis3 > lppt + 3: # 3 means the AG + one nucleotide distance, no U2AF can bind the downstream sequence of BPS if this is true
#            pptS = seq[ipos+lmotif:sL-3] # bps sequence
#            pptSC = bppt_get_pptsc(pptS,lppt,l_max_ppt,baseppt)
#        
#        SC = (bpSC * pptSC)
#        inp = bpS+"\t"+str(pAG)+"\t"+str(bpSC)+"\t"+str(pptSC)+"\t"+str(SC)
#        #print pAG
#       
#        if len(orinp) == 0:
#            orinp.append(inp)
#        else:
#            flag_in = 0
#            for i in range(0,len(orinp)):
#                line = orinp[i].split("\t")
#                scold = float(line[4])
#                if scold < SC:
#                    orinp.insert(i,inp)
#                    flag_in = 1
#                    break
#            if flag_in == 0:
#                orinp.append(inp)
#
#        totsc += SC
#        totbpsc += bpSC
#        totpptsc += pptSC
#        npos += 1
#
#    msc = totsc/npos
#    mbpsc = totbpsc/npos
#    mpptsc = totpptsc/npos
#    dsc = []
#    dbpsc = []
#    dpptsc = []
#    sdsc = 0
#    sdbpsc = 0
#    sdpptsc = 0
#    for i in range(0,len(orinp)):
#        line = orinp[i].split("\t")
#        sc = float(line[4])
#        bpsc = float(line[2])
#        pptsc = float(line[3])
#        dd = sc-msc
#        dsc.append(dd)
#        sdsc += dd*dd
#        dd = bpsc-mbpsc
#        dbpsc.append(dd)
#        sdbpsc += dd*dd
#        dd = pptsc-mpptsc
#        dpptsc.append(dd)
#        sdpptsc += dd*dd
#
#    sdsc = math.sqrt(sdsc/npos)
#    sdbpsc = math.sqrt(sdbpsc/npos)
#    sdpptsc = math.sqrt(sdpptsc/npos)
#    zsc = []
#    zbps = []
#    zppt = []
#    for i in range(0,len(dsc)):
#        zsc.append(dsc[i]/sdsc)
#        zbps.append(dbpsc[i]/sdbpsc)
#        zppt.append(dpptsc[i]/sdpptsc)
#
#    return(orinp,zbps,zppt,zsc)
#
















def bppt_get_pwm(pwmf): # get the pwm of BPS information
    index = 0
    PWMBP = {}
    with open(pwmf,'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            PWMBP[index] = {}
            line = tmp.split('\t')
            PWMBP[index]['A'] = float(line[1])
            PWMBP[index]['C'] = float(line[2])
            PWMBP[index]['G'] = float(line[3])
            PWMBP[index]['T'] = float(line[4])
            index += 1
    return(PWMBP)

def bppt_get_ppt(pptf): # get the ppt score
    PPTS = {}
    with open(pptf,'r') as IN:
        for tmp in IN:
            tmp = tmp.strip()
            if tmp.startswith("#"):
                continue
            line = tmp.split('\t')
            PPTS[line[0]] = float(line[4])
    return(PPTS)
  
def bppt_bpscore(cbps):
    cbpsc = 1
    cbpsN = list(cbps)
    for i in range(0,len(cbpsN)):
        if i in PWMBP:
            cbpsc = cbpsc*PWMBP[i][cbpsN[i]]
        else:
            print("Error: the input bps is longer than the PWM")
    return(cbpsc)

def bppt_get_bpscore(l,basebp):
    nn = ['A','C','G','T'] # nucleotides
    NNS = ['A','C','G','T'] #  nucleotides combination with fixed length
    NN = [] #  nucleotides combination with fixed length

    for i in range(2,l+1):
        for nns in NNS:
            for n in nn:
                newN = nns+n
                NN.append(newN)
        NNS = NN 
        NN = []

    basebpsc = bppt_bpscore(basebp)
    cBPSC = {}
    for nn in NNS:
        cBPSC[nn] = bppt_bpscore(nn)/basebpsc
    
    return(cBPSC)

def bppt_get_AGEZ(seq,offset=12,maxL=-1): # get the AGEZ 
    sL = len(seq)
    ss = seq[0:-offset].split("AG")
    if maxL > 0:
        pAG = sL - maxL
    else:
        pAG = sL - (offset+len(ss[-1])+14)
    if pAG < 0:
        pAG = 0
    return(pAG)

def bppt_get_pptsc(pptS,lppt,l_max_ppt):
    end = len(pptS)-lppt
    if end > l_max_ppt: # if the ppt sequence is larger than l_max_ppt, only check the l_max_ppt
        end = l_max_ppt
    pptbasesc = (l_max_ppt-lppt+1)*PPTS["TTTTTTTT"]
    pptSC = 0
    for i in range(0,end):
        cppts = pptS[i:i+lppt]
        pptSC += PPTS[cppts]
    if pptbasesc == 0:
        return(pptSC)
    else:
        return((pptSC/pptbasesc))

def bppt_dis_pro(pAG,offset=22,interval=4):#4 best
    if abs(pAG-offset) > 700:
        return(0)
    else:
        return(1/(math.exp((abs(pAG-offset))/interval)+1))#best

def bppt_get_BPPTsc(seq,maxL): # get the candidate bps and ppt and their score
    print(maxL, "...", seq)

    baseppt = "TTTTTTTT"
    lmotif = 7
    lppt = 8
    l_max_ppt = 20
    pstart = bppt_get_AGEZ(seq=seq,offset=12,maxL=maxL)
    sL = len(seq)

    totbpsc = 0
    totpptsc = 0
    totsc = 0
    npos = 0

    orinp = []
    orsc = []

    for ipos in range(pstart,sL-3-lmotif): # the BPS could close to the 3' end with only one nucleotide distance 
    #for ipos in range(pstart,sL-14-lmotif): # the BPS could close to the 3' end with only one nucleotide distance 
        pAG = sL - ipos - 5
        bpS = seq[ipos:ipos+lmotif] # bps sequence
        bpSC = cBPSC[bpS] # bps sequence score

        pptSC = 0
        dis3 = pAG - 1 # the distance of BPS last nucleotide to the 3' end
        if dis3 > lppt + 3: # 3 means the AG + one nucleotide distance, no U2AF can bind the downstream sequence of BPS if this is true
            pptS = seq[ipos+lmotif:sL-3] # bps sequence
            pptSC = bppt_get_pptsc(pptS,lppt,l_max_ppt)
        
        SC = (bpSC * pptSC)
        inp = bpS+"\t"+str(pAG)+"\t"+str(bpSC)+"\t"+str(pptSC)+"\t"+str(SC)
        #print pAG
       
        if len(orinp) == 0:
            orinp.append(inp)
        else:
            flag_in = 0
            for i in range(0,len(orinp)):
                line = orinp[i].split("\t")
                scold = float(line[4])
                if scold < SC:
                    orinp.insert(i,inp)
                    flag_in = 1
                    break
            if flag_in == 0:
                orinp.append(inp)

        totsc += SC
        totbpsc += bpSC
        totpptsc += pptSC
        npos += 1

    msc = totsc/npos
    mbpsc = totbpsc/npos
    mpptsc = totpptsc/npos
    dsc = []
    dbpsc = []
    dpptsc = []
    sdsc = 0
    sdbpsc = 0
    sdpptsc = 0
    for i in range(0,len(orinp)):
        line = orinp[i].split("\t")
        sc = float(line[4])
        bpsc = float(line[2])
        pptsc = float(line[3])
        dd = sc-msc
        dsc.append(dd)
        sdsc += dd*dd
        dd = bpsc-mbpsc
        dbpsc.append(dd)
        sdbpsc += dd*dd
        dd = pptsc-mpptsc
        dpptsc.append(dd)
        sdpptsc += dd*dd

    sdsc = math.sqrt(sdsc/npos)
    sdbpsc = math.sqrt(sdbpsc/npos)
    sdpptsc = math.sqrt(sdpptsc/npos)
    zsc = []
    zbps = []
    zppt = []
    if sdsc != 0:
        zsc = [i/sdsc for i in dsc]
    if sdbpsc != 0:
        zbps = [i/sdbpsc for i in dbpsc[:len(dsc)]]
    if sdpptsc != 0:
        zppt = [i/sdpptsc for i in dpptsc[:len(dsc)]]

    return(orinp,zbps,zppt,zsc)

#global variables
REPORTN = 1
bpp_path = Path("bin/BPP")
PWMBP = bppt_get_pwm(bpp_path / "pwmBP_human.txt")
PPTS=bppt_get_ppt(bpp_path / "scPPT_human.txt")
cBPSC = bppt_get_bpscore(7,"TACTAAC")