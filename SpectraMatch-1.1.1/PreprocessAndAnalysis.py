#!/usr/bin/env python
#coding: utf-8 
import sys, optparse, os, tempfile, subprocess, shutil, re, math
import datetime
import pandas as pd
import numpy as np

H_PLUS = 1.007276466#1.007825

'''
Mono-isotopic mass, http://61.50.134.137/mascot/help/aa_help.html
'''
dict_aa_mass = {'A':71.037, 'B':114.535, 'C':103.009, 'D':115.027, 'E':129.043, 'F':147.068,
                    'G':57.021, 'H':137.059, 'I':113.084, 'J':113.084, 'K':128.095, 'L':113.084,
                    'M':131.040, 'N':114.043, 'O':237.148, 'P':97.053, 'Q':128.059, 'R':156.101,
                    'S':87.032, 'T':101.048, 'U':150.954, 'V':99.068, 'W':186.079, 'X':111.0,
                    'Y':163.063, 'Z':128.551}

def stop_err(msg, ret=1):
    sys.stderr.write(msg)
    sys.exit(ret)

def time_now():
    return datetime.datetime.now().strftime('%X') 

def timer(t2, t1):
    return float((t2 - t1).seconds) + (t2 - t1).microseconds / 1000000.0


def getPepFromAllFile():
    #===============================================================================
    # Get Real or Decoy psms  
    #===============================================================================
    flist = []
    with open("./f_list.txt",'r') as f:
        flist = [ x.strip() for x in f ]

    
    foutTarget = './pep_pro_pair_target.txt'
    foutDecoy  = './pep_pro_pair_decoy.txt'
        
    foutT = open(foutTarget,'w')
    foutD = open(foutDecoy,'w')
    
    for index in range( len(flist) ) :#len(flist)

        fullpath = os.path.join('/media/_EXTEND2/cAOrUIfANG1989/percolator_csv/human/QC_input/',flist[index])
        print fullpath
        try:
            f = open(fullpath,'r')
        except:
            continue
        f.readline()
        lineNum = 0
        for line in f:
            lineNum += 1
            line = line.split(',')
            #pep_d = {}
            #pep_d['sequence'] = line[0]
            sequence = line[4]
            proteins = line[5].strip()
            leadings = line[5].strip()
            try:
                #score    = line[46]
                score    = line[2] # PEP
            except:
                continue
            
            if not score or score=='NaN':
                score = '0'
                continue # can't judge
            #pep_d['proteins'] = line[8]
            
            record = str(index)+'|'+str(lineNum)
            
            if 'REV_' in proteins:
                #if ';' in leadings:
                    #print '#' + leadings
                tmp = [ record, sequence, proteins,score]
                foutD.write( '\t'.join(tmp) + '\n')
            elif 'CON_' in proteins or 'CON_' in leadings:
                #if ';' in proteins:
                    #print '#' + leadings
                tmp = [ record, sequence, proteins,score]
                foutD.write( '\t'.join(tmp) + '\n')
            else:
                tmp = [ record, sequence, proteins,score]
                foutT.write( '\t'.join(tmp) + '\n')
            #pep_list.append(pep_d)
        f.close()
        
        
    foutT.close()
    foutD.close()
    #print 'peptide_list.append...[%s]' % type, time_now()
    exit(0)
    #return pep_list

def get_ibaq_num(pro_seq):
    length=len(pro_seq)
    tmp_len=0
    ans=0
    for i in range(length):
        tmp_len=tmp_len+1
        if pro_seq[i]=='K' or pro_seq[i]=='R':
            if i==length-1 or pro_seq[i+1]!='P':
                if tmp_len>=7 and tmp_len<=40:
                    ans=ans+1
            tmp_len=0
    
    if ans == 0 : ans = 1
    return ans      


def get_pro_Length(peptide_list, qvalueThres, ionThres, pro_len_existed):
    dict_pro_length = {}
    tmp_set_pro = set()
    for line in f:
        acc = line[0]
        len = line[1].strip()
        dict_pro_length[acc] = int(len)

    return dict_pro_length

def pep_pro_association(peptide_list, qvalueThres, ionThres):
    tmp_set_pep = set()
    tmp_set_pro = set()
    dict_pep_has_pro = {}
    dict_pro_has_pep = {}
    #pep_index = -1
    nonredundant_pep = 0
    c = 0
    for peptide in peptide_list:
        pep_index = str(c)+'|0'
#         type = peptide['type']
        
        sequence = peptide[0]
        pga_list = peptide[1]
#         if len(sequence) < 7:
#             # print seq
#             continue
        try:
            score = 0
            #print score
        except:
            #print peptide[3]
            score = 0.0
        if score > qvalueThres:
            continue
#         if ionscore < ionThres:
#             continue
#         qvalue = peptide['q_value']
        #=======================================================================
        # Round(ionscore)
        #=======================================================================
        # ionscore = round(ionscore)
#         if qvalue > qvalueThres:
#             continue
        #nonredundant_pep += 1
        #area = peptide['area']
        
        
        #charge = peptide['charge']
        # seq = seq + '_' + str(ionscore)
        # seq = seq + '_' + str(charge)
        seq = sequence.upper()
        if seq in dict_pep_has_pro:
            #===============================================================
            """ Maybe many peptides with same ionscore """
            #===============================================================
            #dict_pep_has_pro[seq]['ionscore'] += ionscore
            dict_pep_has_pro[seq][2] += score # score
            #dict_pep_has_pro[seq]['area'] += area
            dict_pep_has_pro[seq][1] += 1 # psms
            #pass
            
        else:
            #tmp_set_pep.add(seq)
            tmp = [set(), 1, score] # proteins,psms,score
            
            #tmp['ionscore'] = ionscore
            #tmp['charge'] = charge
            #tmp['area'] = area
            
            dict_pep_has_pro[seq] = tmp
        #=======================================================================
        """ Get every protein in which this peptide appears """
        #=======================================================================
        
        for pro in pga_list.split(';'):

            dict_pep_has_pro[seq][0].add(pro)
            
            if pro in dict_pro_has_pep:     
                #dict_pro_has_pep[pro]['peptides'].add(seq)
                dict_pro_has_pep[pro][0].add(seq)
                #dict_pro_has_pep[pro]['pep_index'].add(pep_index)
                dict_pro_has_pep[pro][1].add(pep_index)
                
            else:   
                symbol = ''  # getSymbol(pro)
                #tmp_set_pro.add(pro)
                tmp = [set([seq]), set([pep_index]), set() ] # peptides, pep_index, sameset
                #tmp['peptides'] = set([seq])
                #tmp['pep_index'] = set([pep_index])
                #tmp['sameset'] = set()
                #tmp['symbol'] = symbol
                
                #tmp['type'] = type
                
                dict_pro_has_pep[pro] = tmp
    for x in dict_pep_has_pro:
        s = dict_pep_has_pro[x][2]
        psms = dict_pep_has_pro[x][1]
        dict_pep_has_pro[x][2] = s/psms
        
    return (dict_pep_has_pro, dict_pro_has_pep, nonredundant_pep)

def getProLength(pro):
    try:
        L = int( file('./proLength/%s'%pro).read().strip() )
    except:
        #print pro,":length err"
        L = 100000
        
    return L

def useParsimony(eid,type, qvalueThres, ionThres, peptide_list, flag, dict_pro_length):
    flag = 1
    t1 = datetime.datetime.now()
    if flag != 0:
        print 'Grouping...[%s]' % type, t1.strftime('%X')
    
    (dict_pep_has_pro, dict_pro_has_pep,nonredundant_pep) = pep_pro_association(peptide_list, qvalueThres, ionThres)
    
    #===========================================================================
    # print 'for x,y in dict_pep_has_pro.iteritems()'
    # i=0
    # for x,y in dict_pep_has_pro.iteritems():
    #     i+=1
    #     if i>10:break
    #     print x,y
    # print 'for x,y in dict_pro_has_pep.iteritems()'    
    # i=0    
    # for x,y in dict_pro_has_pep.iteritems():
    #     i+=1
    #     if i>10:break
    #     print x,y
    #     
    # print 
    # exit(0)
    #===========================================================================
    # f = open('outpep%s(%s).txt' % (eid, type), 'w')           
    set_uni_pep = set()
    for seq, info in dict_pep_has_pro.iteritems():
        # f.write(seq.split('_')[0] + '\n')
        # f.write(seq.split('_')[0].upper() + '\t'+';'.join(list(set(list_pro)))+'\n')
        if len(info[0]) == 1:# info['proteins'] is a set
            set_uni_pep.add(seq)
    # f.close()
    if flag != 0:
        print 'len(dict_pep_has_pro):', len(dict_pep_has_pro)
        print 'len(dict_pro_has_pep):', len(dict_pro_has_pep)
        print 'len(set_uni_pep):', len(set_uni_pep)

 
    #===============================================================================
    ''' Get distinct proteins identified by only discrete peptide(s)'''
    #===============================================================================
    j = k = 0
    dict_uni_pro      = {} # Proteins have unique peptides
    dict_none_uni_pro = {} # Proteins don't have unique peptides
    pro_onlyUniPep    = set()
    for pro, info in dict_pro_has_pep.iteritems():
        proscore = 0
        pep_psms = 0
        peptides = info[0]  # Is a set
        pepnum = len(peptides)
        pep_index_set = info[1] # Is a set
        symbol = ''#info['symbol']
        # peps_other = list(peptides)
        uni_peps = []
        has_uni = 0
        # print pepnum

        for peptide in peptides:
            ''' Pro_score comes from ionscore '''
            proscore += dict_pep_has_pro[peptide][2]
            pep_psms += dict_pep_has_pro[peptide][1]
            if peptide in set_uni_pep:
                has_uni += 1
                uni_peps.append(peptide)
                # peps_other.remove(peptide)
        sameset = set()
        tmp = [ peptides, pep_index_set, sameset, uni_peps, symbol, proscore, pep_psms ]
#         tmp['peptides'] = peptides  # Is a set
#         tmp['pep_index'] = pep_index_set # Is a set
#         tmp['sameset'] = set()
#         tmp['uni'] = uni_peps
#         tmp['symbol'] = symbol
#         tmp['score'] = proscore
#         tmp['pep_psms'] = pep_psms
#         tmp['type'] = type
#         tmp['score'] = proscore
#         tmp['pep_psms'] = pep_psms
        if has_uni > 0:
            dict_uni_pro[pro] = tmp
            j += 1
            if has_uni == pepnum:
                pro_onlyUniPep.add(pro)
        else:
            #tmp['score'] = 0
            dict_none_uni_pro[pro] = tmp
            
    if flag != 0:
        print 'Finished getting dict_uni_pro...[%s]' % type, time_now()
        print 'dict_none_uni_pro:', str(len(dict_none_uni_pro))
    # exit(0)
    
    ''' Start to get sameset, and merge to dict_uni_pro '''
    sameset_find = total_sameset = 0
    searched = set()
    pro_sameset = set()
    dict_sameset_pro = {}
    for pro1, info1 in dict_none_uni_pro.items():
        if pro1 in searched : continue
        searched.add(pro1)
        master_pro = pro1
        #len_pro1 = get_pro_len_Simple(conn, meta, pro1)
        len_pro1 = getProLength(pro1)
        min_len = len_pro1
        sameset_find = 0
        tmp = []
        for pro2, info2 in dict_none_uni_pro.items():
            if pro2 in searched : continue
            if info1[0] == info2[0]:  # info[0] == set(peptides) #Is a set
                #len_pro2 = get_pro_len_Simple(conn, meta, pro2)
                len_pro2 = getProLength(pro2)
                if len_pro2 < min_len:
                    min_len, master_pro = len_pro2, pro2
                    
                symbol = ''#info2['symbol']
                sameset_find += 1
                # print pro1, pro2
                if sameset_find == 1: 
                    del dict_none_uni_pro[pro1]
                    del dict_none_uni_pro[pro2]
                    searched.add(pro2)
                    #     [ peptides, pep_index_set, sameset,           uni_peps, symbol, proscore,  pep_psms ]
                    tmp = [ info1[0], info1[1],      set([pro1, pro2]), [],       symbol, info1[5],  info1[6] ]
#                     tmp['peptides'] = info1['peptides']
#                     tmp['pep_index'] = info1['pep_index']
#                     tmp['sameset'] = set([pro1, pro2])
#                     tmp['uni'] = []
#                     tmp['symbol'] = symbol
#                     tmp['score'] = info1['score']
#                     tmp['pep_psms'] = info1['pep_psms']
#                     tmp['modification'] = info1['modification']
#                     tmp['type'] = type
                else:
                    # pro_sameset.add(pro2)
                    del dict_none_uni_pro[pro2]
                    searched.add(pro2)
                    tmp[2].add(pro2) # sameset
                    tmp[1] = tmp[1] | info2[1]# pep_index
                # dict_uni_pro[pro_uni]['sameset'].append(pro)
        # print sameset_find,'====================='
        if sameset_find > 0:
            total_sameset += 1
            tmp[2].remove(master_pro)
            dict_sameset_pro[master_pro] = tmp
    
    def isUniq(pro,pep):
        for key,info in dict_uni_pro.iteritems():
            if key == pro:continue
            if pep in info[0]:
                return False
            
        for key,info in dict_sameset_pro.iteritems():
            if key == pro:continue
            if pep in info[0]:
                return False
            
        for key,info in dict_none_uni_pro.iteritems():
            if key == pro:continue
            if pep in info[0]:
                return False
            
        return True
    for key,info in dict_sameset_pro.items():
        setPeptides = info[0]
        for p in setPeptides:
            if isUniq(key,p):
                info[3].append(p)
            
    dict_uni_pro.update(dict_sameset_pro)
    
    """ Start to get remaining protein, and merge to dict_uni_pro """
    remaining = 0
    #print 'Remaining:\n'
    for pro, info in dict_none_uni_pro.items():
        equal_cover = 1
        pep_num = len(info[0])
        for pep in info[0]:
            if equal_cover == 0:break
            for next_pro in dict_pep_has_pro[pep][0]:
                if equal_cover == 0:break
                if next_pro == pro:continue
                pep_num_next_pro = len(dict_pro_has_pep[next_pro][0])
                if pep_num < pep_num_next_pro :
                    equal_cover = 0
        if equal_cover:
            #print pro
            remaining += 1
            dict_uni_pro[pro] = info
            del dict_none_uni_pro[pro]
              
    if flag != 0:
        print 'Pro has unipep:', j
        print 'Pro has Only unipep:', len(pro_onlyUniPep)
        print 'Pro has sameset:', total_sameset
        print 'Assigned Pro:', len(dict_uni_pro)
        print 'Remaining Pro:',remaining
        print 'Unassigned:', len(dict_none_uni_pro)-remaining
        print 'Finish useParsimony(%s,  %s)..................................' % (eid, type), time_now()
    
    if type == 1:
        return (dict_uni_pro, dict_none_uni_pro, dict_pep_has_pro, dict_pro_has_pep, nonredundant_pep)
    else:
        return (dict_uni_pro, nonredundant_pep)

    
def startProAssemble(eid_list, peptide_list, fdrThres, ionThres):
    f_result = open('./result.txt','a')
    #dir = '/usr/local/galaxyDATA01/incubator/python/ms_parser'
    eid = eid_list[0]
    t1 = datetime.datetime.now()
    print '\nStart startProteinAssemble...eid =', eid
    f_result.write('start'+'\n')
    f_result.flush()
    #ionThres = 20
    step = 0.0001
    step_score = 1
    step_pep = 0.0001
	
    fdrThres = fdrThres + 0.0005
	
    pro_score_Thres = 20
    qvalueThres = 0.001
    pepThres    = 0.0032
	
    print '\n# # # # # Ionscore Threshold:', ionThres
    #===========================================================================
    # Start cycle to find FDR < 0.01
    #===========================================================================
    FDR = FDR_pep = FDR_psms = 1
    times = 0
    flag = 0
    dict_pro_length = {}
    (dict_uni_pro_real, dict_none_uni_pro_real, dict_pep_has_pro_real, dict_pro_has_pep_real, total_pep_1) = useParsimony(eid, 1, pepThres, ionThres, peptide_list, flag, dict_pro_length)
    psms_fdr_real = psms_fdr_decoy = 0
    del peptide_list
    t2 = datetime.datetime.now()
    delta = timer(t2, t1)  # float((t2 - t1).seconds) + (t2 - t1).microseconds / 1000000.0
    print '\nGroup Cycle used %.3f seconds' % delta
    #===============================================================================
    # Output proteins
    #===============================================================================
    print 'Protein group_result...' , time_now()
    
    #exit(0)
    
    
    tmp_set_pro = set()
    tmp_set_pep_index = set()
    tmp_set_pep_pga = set()
    set_pep_index = {}
    pga_peptide = {}
    group_result = []   
    total_FOT = 0.0
    for pro, info in dict_uni_pro_real.iteritems():
        #=======================================================================
        """ Using Pep_index to know which peptide used to assemble proteins """
        #=======================================================================
        pep_index_set = info[1]
        tmp_set_pep_index = tmp_set_pep_index | pep_index_set

        #=======================================================================
        """ Below is getting Protein info """
        #=======================================================================
        proscore = info[5]
        symbol = info[4]
        pep_psms = info[6]
        psms_fdr_real += pep_psms
        
        sameset = ';'.join(sorted(info[2]))
        peps = sorted(info[0])
        num_pep = len(peps)
        uni_peps = sorted(info[3])
        num_unipep = len(uni_peps)
        area_pro = 0.0
        #=======================================================================
        """ Distribute peptide area to Protein """
        #=======================================================================
        for peptide in peps:
            if peptide in tmp_set_pep_pga:
                pga_peptide[peptide].add(pro)
            else:
                tmp_set_pep_pga.add(peptide)
                pga_peptide[peptide] = set([pro])
#                 
#             assembled = 0
#             area_pep = dict_pep_has_pro_real[peptide]['area']
#             pros_of_pep = dict_pep_has_pro_real[peptide]['proteins']
#             for p in pros_of_pep:
#                 if p in dict_uni_pro_real: assembled += 1
#             area_average = area_pep / assembled if assembled != 0 else 0
#             area_pro += area_average
            
        tmp = [pro, symbol, proscore, sameset, num_pep, num_unipep, pep_psms, area_pro, ';'.join(sorted(pep_index_set))]
#         tmp['accession'] = pro
#         tmp['symbol'] = symbol
#         tmp['score'] = proscore
#         tmp['sameset'] = sameset
#         tmp['num_pep'] = num_pep
#         tmp['num_unipep'] = num_unipep
#         tmp['pep_psms'] = pep_psms
#         tmp['area'] = area_pro
#         tmp['pep_index'] = pep_index_set
        #tmp['modification'] = modification
        group_result.append(tmp)
        
    #===========================================================================
    # Print proteins with proper FDR        
    #===========================================================================
    
    output_nums = 0
    group_result.sort(key=lambda x:x[2], reverse=1)
 
    output_nums = len(group_result)
    #===============================================================================
    # Unassigned proteins
    #===============================================================================
    print 'Unassigned proteins:',len(dict_none_uni_pro_real)
    
#    FDR_psms = float(psms_fdr_decoy)/psms_fdr_real if psms_fdr_real !=0 else 1
    
#    print 'FDR of psms =',FDR_psms
#    print 'FDR of peptide =',FDR_pep
    print '\nFinal pep Threshold = ', pepThres
#    pepfinal = open('./pepfinal','w')
#    pepfinal.write(str(pepThres))
#    pepfinal.close()
    print 'After cycled %s times\nWrite %s proteins\nWith FDR = %f and ionscoreThreshold = %s......' % (times, output_nums, FDR, ionThres), time_now()
    #===============================================================================
    # Time used
    #===============================================================================
    t2 = datetime.datetime.now()
    delta = timer(t2, t1)  # float((t2 - t1).seconds) + (t2 - t1).microseconds / 1000000.0
    print 'Over...', str(t2.strftime('%X'))
    print 'Used %.3f seconds' % delta
    print 'EID:', eid
    print 
    
    writePepIndex(tmp_set_pep_index)
    print 'writePepIndex()  done'
    
    writePga_peptide(pga_peptide)
    print 'writePga_peptide()  done'
    
    writeProtein(group_result)
    print 'writeProtein() done'
    f_result.close()
    return (group_result, qvalueThres, FDR, dict_uni_pro_real, tmp_set_pep_index, pga_peptide, dict_pro_length)

def getPeps(type,file_peptide):
    #foutTarget = file_peptide #'./pep_pro_1.txt'
    foutDecoy  = './pep_pro_pair_decoy.txt'
    pepList = []
    if type == 1:
        f=open(file_peptide,'r')
        for line in f:
            line = line.split('\t')
#            fileIndex = line[0]
            seq = line[0]
            pros = line[1].strip()
#            score = line[3].strip()
            pepList.append([seq,pros])
        f.close()
            
    else:
        f=open(foutDecoy,'r')
        for line in f:
            line = line.split('\t')
            fileIndex = line[0]
            seq = line[1]
            pros = line[2].strip()
            score = line[3].strip()
            if ';' in pros:
                tmp = []
                #print pros
                for p in pros.split(';'):
                    if 'REV_' in p:
                        tmp.append(p.split('REV_')[1])
                    elif 'CON_' in p:
                        tmp.append(p.split('CON__')[1])
                    else:
                        tmp.append(p)
                pros = ';'.join(tmp)
            else:
                if 'REV_' in pros:
                    pros = pros.split('REV_')[1]
                elif 'CON_' in pros:
                    pros = pros.split('CON_')[1]
                    
                
            pepList.append([fileIndex,seq,pros,score])
        f.close()
    
    return pepList
    
def cacheProtein( eid_list, file_data, stype, fdrThres, ionThres, file_peptide):
    # updateExpStage(conn,meta,file_data)
    #t1 = datetime.datetime.now()
    print '\nStart cacheProtein...ionThres=',ionThres

    rid = file_data['r_num']
    rank = file_data['rank']
    #print '[Exp ID %s] [Rep %s] [Rank %s] [Type %s]' % (eid_list[0], rid, rank, stype)
    
    sid = 0
    
    peptide_list   = getPeps(1,file_peptide)
    #pep_list_decoy = getPeps(-1)

    #peptide_list = []        
    (group_result, qvalueThres, fdr, dict_uni_pro, set_pep_index, pga_peptide, dict_pro_length) = startProAssemble(eid_list, peptide_list, fdrThres, ionThres)
   
    
        # file_data, stype, qvalueThres, , , redundant_pep, set_pep_index, pga_peptide
    return ( qvalueThres, set_pep_index, pga_peptide)

def writeProtein(group_result):
    f = open('pep_qc_proteinAll.txt','w')
    # p = [pro, symbol, proscore, sameset, num_pep, num_unipep, pep_psms, area_pro, pep_index_set]
    for p in group_result:
        tmp = [ str(x) for x in p ]
        f.write('\t'.join(tmp) + '\n')
    f.close()
    
def writePga_peptide(pga_peptide):
    f = open('pep_qc_pga_peptide.txt','w')
    for i,j in pga_peptide.iteritems():# set_pep_index is a set() !!!!!!
        #fileNumber, lineNumber = i.split('|')
        f.write(i + '\t' + ';'.join(sorted(j)) + '\n')
        #pep2print[ int(fileNumber) ].append(int(lineNumber))
    f.close()
    
def writePepIndex(set_pep_index):
    #
    f = open('pep_qc_pep2print.txt','w')
    for i in set_pep_index:# set_pep_index is a set() !!!!!!
        #fileNumber, lineNumber = i.split('|')
        f.write(i + '\n')
        #pep2print[ int(fileNumber) ].append(int(lineNumber))
    f.close()
    
def cachePeptide(file_data, stype, qvalueThres, ionThres, set_pep_index, pga_peptide) :
    
    #redundant_pep = getPepFromDB(eid, rid, rank, 1, stype, ionThres)
    a = []
    total_FOT = 0.0
    pep_has_pro = {}
    set_tmp = set()
    pep2print = [ [] for i in range(78) ]
    f = open('pep_qc_pep2print.txt','r')
    for i in set_pep_index:# set_pep_index is a set() !!!!!!
        fileNumber, lineNumber = i.split('|')
        pep2print[ int(fileNumber) ].append(int(lineNumber))
    f.close()
        
def my_QC(options, file_peptide):

    print 'Start caching Exp...', time_now()    
    #exit(0)
    
    file_data = {}
    file_data['r_num'] = 1
    file_data['rank'] = 1
    
    fdrThres = options['min_fdr']
    ionThres = options['min_ion']
    eid_list = options['exp_list']
    
    stype = options['stype']
    
    ( qvalueThres, set_pep_index, pga_peptide) = cacheProtein(eid_list, file_data, stype, fdrThres, ionThres, file_peptide)
    
    peptides = cachePeptide(file_data, stype, qvalueThres, ionThres, set_pep_index, pga_peptide)      
     
    print 'Complete caching Exp...', time_now()
    
    #return (proteins, peptides, genes)
	
#ConvrtFile2Txt
def GetNum(String):
    return int(re.findall(r'(\w*[0-9]+)\w*',String)[0])
# =============================================================================
def MethodSelect(TxtFilePath,mtStr,RootFilePath,SortedList):
	TxtFileDenoisePath=TxtFilePath+mtStr
	os.mkdir(RootFilePath+TxtFileDenoisePath)
	TxtDenoiseName=TxtFileDenoisePath+'.txt'
	DenoiseFilePath=RootFilePath+TxtFileDenoisePath+"/"
	f_TxtDenoiseName=open(TxtDenoiseName,'w')
	for i in range(0,len(SortedList)):
		f_TxtDenoiseName.write(DenoiseFilePath+"scan="+str(SortedList[i])+".txt"+'\n')
	f_TxtDenoiseName.close()
# =============================================================================
def ConvertRawFile2Txt(InputFile,OutputFile,method='all'):
	# ==============================================
	RootFilePath="./"
	TxtFilePath=OutputFile
	os.mkdir(TxtFilePath)
	NameBeforeDenoise=OutputFile+'.txt'
	# ==============================================
	with open(InputFile,'r') as f_InputFile:
		DataFromInputFile=f_InputFile.read()
	DataSplit=DataFromInputFile.split("\nEND IONS\n")
	NameList=[]
	for Datasplit in DataSplit[:-1]:
		GetLines=Datasplit.split("\n")
		if GetLines[1].startswith('TITLE'):
			NameGet = GetLines[1].split(' ')[-1].replace('''"''',"")
		NameList.append(NameGet)
		PeptideMass=GetLines[3].replace("PEPMASS=","").strip()
		PeptideMass=PeptideMass.split(' ')
		JudgeCharge=0
		if (len(GetLines)>4):
			if GetLines[4].startswith('CHARGE'):
				PeptideCharge=GetLines[4].split('=')[1].split('+')[0]
				JudgeCharge=5
			else:
				PeptideCharge=1
				JudgeCharge=4
			with open(TxtFilePath+"/"+NameGet+".txt",'w') as f1_NameGet:
				for Getlines in GetLines[JudgeCharge:]:
					PosValue=Getlines.split(" ")
					f1_NameGet.write(str(PosValue[0])+'\t'+PosValue[1]+'\r\n')
				f1_NameGet.write(str(PeptideMass[0])+'\t'+str(PeptideCharge))
		else:
			continue
	# ==============================================
	TxtOutFilePath=RootFilePath+TxtFilePath+"/"
	TxtFileList=os.listdir(TxtFilePath+"/")
	# ==============================================
	NumList=[]
	[NumList.append(GetNum(TxtFileList[i])) for i in range(0,len(TxtFileList))]
	SortedList=sorted(NumList)
	f_NameBeforeDenoise=open(NameBeforeDenoise,'w')
	for i in range(0,len(SortedList)):
		f_NameBeforeDenoise.write(TxtOutFilePath+"scan="+str(SortedList[i])+".txt"+'\n')
	f_NameBeforeDenoise.close()
	# ==============================================
	if method=='CutThreshold':
		MethodSelect(TxtFilePath,'_CutThreshold',RootFilePath,SortedList)
	if method=='ZoneSelect':
		MethodSelect(TxtFilePath,'_ZoneSelect',RootFilePath,SortedList)
	if method=='ChargeSelect':
		MethodSelect(TxtFilePath,'_ChargeSelect',RootFilePath,SortedList)
	if method=='all':
		MethodSelect(TxtFilePath,'_CutThreshold',RootFilePath,SortedList)
		MethodSelect(TxtFilePath,'_ZoneSelect',RootFilePath,SortedList)
		MethodSelect(TxtFilePath,'_ChargeSelect',RootFilePath,SortedList)

#sort_output_file
def sort_output(file1,file2,file3,file4,file5,thres=5):
    os.mkdir(file4) #会有输出
    df=pd.read_csv(file1,sep="\t",header=None,index_col=None)
    #df[8]=list(map(lambda x:1 if x[0].split('/')[-1].split('_')[0]==x[1].split('/')[-1].split('_')[0] else 0,df[[0,1]].values))
    df[8]=list(map(lambda x:float(x[0])+float(x[1])+2*float(x[2]),df[[2,3,4]].values))
    frames=[]
    sample=[]
    length_df=len(df.index)
    length_sa=len(set(df[0]))
    length_df_1=len(df.index)-1
    start=0
    end=0
    cnt=0
    for i in range(length_sa):
        frames.append('df'+str(i))
    f=open(file3,'w')
    for i in range(length_df):
        s=df.loc[i,0].split('/')[-1]
        f.write(s+'\n')
        if s not in sample:
            end=i
            if end != 0 and len(sample)==cnt+1:
                #print ("555")
                df1=df.loc[start:(end-1),:]
                df1=df1.sort(columns=[8],ascending=0)
                df1=df1[[0,1,2,3,4,5,6,8,7]]
                df1.to_csv(file4+'/'+str(sample[-1]),sep='\t',header=None,index=None)
                frames[cnt]=df1[:5]
                cnt=cnt+1
                del df1
            start=i
            sample.append(s)
            if end==length_df_1:
                df1=df.loc[end:end,:]
                df1=df1.sort(columns=[8],ascending=0)
                df1=df1[[0,1,2,3,4,5,6,8,7]]
                df1.to_csv(file4+'/'+str(sample[-1]),sep='\t',header=None,index=None)
                frames[cnt]=df1[:5]
                cnt=cnt+1
                del df1
        else:
            end=i
            if end==length_df_1:
                df1=df.loc[start:end,:]
                df1=df1.sort(columns=[8],ascending=0)
                df1=df1[[0,1,2,3,4,5,6,8,7]]
                df1.to_csv(file4+'/'+str(sample[-1]),sep='\t',header=None,index=None)
                frames[cnt]=df1[:5]
                cnt=cnt+1
                del df1
    f.close()
    result=pd.concat(frames)
    result=result[[0,1,2,3,4,5,6,8,7]]
    df_peptide=result[7]
    result.to_csv(file2,sep='\t',header=None,index=None)
    df_peptide.to_csv(file5,header=None,index=None)

#pep_pro_1
def generate_pep_pro(fasta_file,peptide_file,pep_pro_file):
	f_fasta = open(fasta_file,'r')
	f_pep = open(peptide_file,'r')
	pro = ''
	seq = ''
	pro_seq = {}
	for each in f_fasta:
		if '>' in each:
			pro_seq[pro] = seq
			pro = ''
			seq = ''
			pro = each.strip().split('|')[1]
		else:
			seq = seq + each.strip()
	f_out = open(pep_pro_file,'w')
	for each in f_pep:
		pro_set = set()
		for ea in pro_seq:
			if each.strip() in pro_seq[ea]:
				pro_set.add(ea)
		f_out.write(each.strip()+'\t'+';'.join(pro_set)+'\n')
	f_fasta.close()
	f_pep.close()
	f_out.close()

#create_new_lib
def create_new_lib(csv_file,filteredPath,thresHigh=30.0):
	os.mkdir('High-'+filteredPath+'/')
    	path=filteredPath+'/'
    	path1='High-'+filteredPath+'/'
	with open(csv_file,'r') as f:
		data=f.readlines()
	count=0
	for i in data:
		count=count+1
		if 'prot_hit_num' in i:
			break
	dff=pd.read_csv(csv_file,skiprows=count-1,sep=',')
	dff=dff[dff['pep_score'].apply(lambda x:float(x)>=thresHigh)]
	dff2=pd.DataFrame(dff)
	title=dff2.loc[:,'pep_scan_title']
	seq=dff2.loc[:,'pep_seq']
    	title=np.array(title)
    	seq=np.array(seq)
	for i in range(0,len(title)):
		index_file=title[i].split(' ')[-1].replace('''"''',"")
		index_file=index_file.split('~')[0]
		index_file=index_file+'.txt'
		if os.path.exists(path+index_file):
			with open(path+index_file,'r') as f:
				data_got=f.read()
		final_file=seq[i]
		while os.path.exists(path1+final_file+'_'+index_file):
			final_file=final_file+'_'+str(i)+'_'
		with open(path1+final_file+'_'+index_file,'w') as f1:
			f1.write(data_got)

def option_manual():
    return """
Program: PreprocessAndAnalysis.py (Tools for preprocess and analusis in the peptide searching library)
Version: 1.1.1 

Usage:   python PreprocessAndAnalysis.py <command> [options]
Commands options:
        conv
		-raw	  the raw file of lots of spectra in mgf file
		-out      the base folder name for contain more single spectra from mgf file
		-method   CutThreshold,ChargeSelect,ZoneSelect,default All
	sort
		-in       the scoring file from software SpectraMatch
		-out      the selected scoring file by sorting
		-nl       the all peptide name from scoring file
		-fd       the folder to contain every peptide searching result
		-pep      the peptide name to be stitched of selected scoring file by sorting 
		-thres    select scoring file thres,default(5)
	pep_pro
		-fasta    the reference fasta for stitching peptides
		-peptide  the peptides name list to be stitched
		-peppro   the peptides list refering protein
	create
		-csv      the mascot result which contains bing.csv file
        	-filter   the high quality new library refering which filtered folder
		-thres    set which threshold defines high quality
    """
			
def __main__() :
	#convert
	if "conv" in sys.argv:
		InputItem = [item for item in sys.argv if (item.lower()).endswith('.mgf')]
		if "-raw" in sys.argv:
			if len(InputItem) == 1:
				if "-out" in sys.argv:
					Outputfile = sys.argv[sys.argv.index("-out")+1]
					if "-method" in sys.argv:
						method = str(sys.argv[sys.argv.index("-method")+1])
						ConvertRawFile2Txt(InputItem[0],Outputfile,method)
                    			else:
				        	ConvertRawFile2Txt(InputItem[0],Outputfile)
		
	#sort_output_file
	if 'sort' in sys.argv:
		if '-in' in sys.argv:
			index1 =  int((sys.argv).index("-in")) + 1
			source_file1 = sys.argv[index1]
			if '-out' in sys.argv:
				index2 =  int((sys.argv).index("-out")) + 1
				source_file2 = sys.argv[index2]
				if '-nl' in sys.argv:
					index3 =  int((sys.argv).index("-nl")) + 1
					source_file3 = sys.argv[index3]
					if '-fd' in sys.argv:
						index4 =  int((sys.argv).index("-fd")) + 1
						source_file4 = sys.argv[index4]
						if '-pep' in sys.argv:
							index5 =  int((sys.argv).index("-pep")) + 1
							source_file5 = sys.argv[index5]
							if '-thres' in sys.argv:
								index6 =  int((sys.argv).index("-thres")) + 1
								source_file6 = int(sys.argv[index6])
								sort_output(source_file1,source_file2,source_file3,source_file4,source_file5,source_file6)
                            				else:
							     sort_output(source_file1,source_file2,source_file3,source_file4,source_file5)
	
	#pep_pro_1
	if 'pep_pro' in sys.argv:
		if '-fasta' in sys.argv:
			index1 =  int((sys.argv).index("-fasta")) + 1
			source_file1 = sys.argv[index1]
			if '-peptide' in sys.argv:
				index2 =  int((sys.argv).index("-peptide")) + 1
				source_file2 = sys.argv[index2]
				if '-peppro' in sys.argv:
					index3 =  int((sys.argv).index("-peppro")) + 1
					source_file3 = sys.argv[index3]
					generate_pep_pro(source_file1,source_file2,source_file3)
					#getPepFromAllFile() # Run only first time
					
					options = {}
					options['min_fdr'] = 0.01
					options['min_ion'] = 0.01
					options['max_hit'] = 10
					options['dMZ'] = 10
					options['dRT'] = 60
					
					options['exp_list'] = [000000]#[484, 485, 391, 392, 428, 429]
					
					options['stype'] = 'exp'
					my_QC(options, source_file3)
					
					print 'Done'
					#     print peptides[0]
					#     print proteins[0]
					#     print genes[0]
					#     print 'peptides=',len(peptides)
					#     print 'proteins=',len(proteins)
					#     print 'genes=',len(genes)
	#create new library
	if 'create' in sys.argv:
		if '-csv' in sys.argv:
			index1 =  int((sys.argv).index("-csv")) + 1
			source_file1 = str(sys.argv[index1])#sys.argv[index1]
			if '-filter' in sys.argv:
				index2 =  int((sys.argv).index("-filter")) + 1
				source_file2 = str(sys.argv[index2])#sys.argv[index2]
				if '-thres' in sys.argv:
					index3 =  int((sys.argv).index("-thres")) + 1
					source_file3 = float(sys.argv[index3])#sys.argv[index3]
					create_new_lib(source_file1,source_file2,source_file3)
				else:
					create_new_lib(source_file1,source_file2)
	if '-h' in sys.argv:
		print option_manual()
    
if __name__ == "__main__": __main__()


