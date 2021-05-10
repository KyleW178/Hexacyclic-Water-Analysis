# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 23:32:55 2021

@author: Haworth_Lab2
"""

import os,math,csv,pytz,itertools
import numpy as np
from matplotlib import pyplot as PLT
from matplotlib import cm as CM
import seaborn as sns
        
class HP_Interaction():
    def __init__(self):
        self.Carb = ''
        self.object = []
        pass

def aliphatic_C(carbon):
    if carbon[2][0] == 'C' and len(carbon[2]) > 1:
        return True
    else:
        return False 

def water_O(oxygen):
    if oxygen[2][0] == 'O' and oxygen[3] == 'WAT':
        return True
    else:
        return False
#Extract num from atm name
def get_num(atm:str):
    num = ''
    for alp in atm:
        if alp.isdigit():
            num+= alp
    return int(num)

def get_str(atm:str):
    string = ''
    for alp in atm:
        if alp.isalpha():
            string+= alp
    return string


def OO_Clash(c1,c2,crit = 3.5):
    distance = math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)
    if distance < crit:
        return str(round(distance, 3))
    else:
        return False


def HP_HB_Interaction(coo1,coo2,crit,falloff):#For crit =3.5(HP), falloff 0.7. Crit 3.0(HB), falloff 0.5 or other values of the same ratio
    distance = math.sqrt((coo1[0]-coo2[0])**2 + (coo1[1]-coo2[1])**2 + (coo1[2]-coo2[2])**2)
    if (crit - falloff)<= distance <= (crit + falloff):#HP crit 3.5, HB crit 3
        return str(round(distance, 3))
    else:
        return False


#Name: N01, O04, so on
#Consider: regular+ O-O clash(return 0)
def Intra_HB_Crit(name1, name2, aa = 6):
    a, b = get_num(name1), get_num(name2)
    if abs(a-b) in range(2, aa-1):
        if (name1[0]  == 'O' and name2[0] == 'O'):
            return '0'
        else:
            if name1[0] == 'N':
                return '1'
            if name1[0] == 'O':
                return '2'
    else:
        return 0


    
def Intra_HB_Geo(vec1,vec2, limit= 45):
    #VEC1 should be N to H, and VEC2 should be C to O
    unit1 = vec1 / np.linalg.norm(vec1)
    unit2 = vec2 / np.linalg.norm(vec2)
    ang = math.degrees(np.arccos(np.clip(np.dot(unit1, unit2), -1.0, 1.0)))
    final = round(min(ang,(180-ang)),2)

    if  final > limit:
        return 'NG'
    else:
        return final

def retrieve_key(data:dict,compare):
    for keys, values in data.items():
        if values == compare:
            break
    return keys

#    
#def trace_rename_atm(pre_PDB:str,solv_PDB:str):
#    OrigLines = open(pre_PDB,'r').readlines()
#    
#    with open(pre_PDB) as f:
#        for index, l in enumerate(f):
#            pass
#    index = index + 1
#    
#    HB_atm, NH,CO = {} , {}, {} 
#    Alip_C_BB, Alip_C_Exp = {} , {}
#    
#    HB_atm_pair = {}#Gotta figure this one out
#    
#    for atm in OrigLines:
#        split = atm.split()
#        if get_str(split[2]) in ['O' ,'N']: #We don't consider other O or N in this case
#            HB_atm[split[2]] = [split[6],split[7],split[8]]
#            HB_atm_pair[split[2]] = 'NaN'#Indicate which residue 
#        elif split[2][0] == 'C':
#            if get_str(split[2]) == 'C':
#                CO[split[2]] = [split[6],split[7],split[8]]
#            elif get_str(split[2]) == 'CA':
#                Alip_C_BB[split[2]] = [split[6],split[7],split[8]]
#            else:
#                Alip_C_Exp[split[2]] = [split[6],split[7],split[8]]
#        elif split[2][0] == 'H':
#            NH[split[2]] = [split[6],split[7],split[8]]
#    
#    #Update HB_atm pairs
#    for key in HB_atm.keys():
#        if get_str(key) == 'O':#Assume we don't get any hydroxyl or amine right now
#            for Ckey in CO.keys():
#                if get_num(Ckey) == get_num(key):
#                    HB_atm_pair[key] = Ckey
#                    break
#        elif get_str(key) == 'N':
#            for Hkey in NH.keys():
#                if get_num(Hkey) == get_num(key):
#                    HB_atm_pair[key] = Hkey
#                    break
#    
##    print(Alip_C_BB)
##    print(Alip_C_Exp)
##    print(CO)
##    print(HB_atm_pair)
#
#    NewLines = open(solv_PDB,'r').readlines()
#    for n in range(0,index):        
#        newsplit = NewLines[n].split()
#        atmname = newsplit[2]
#        atmcoor = [newsplit[6],newsplit[7],newsplit[8]]
#        if get_str(newsplit[2]) in ['O' ,'N']:
#            
#            HB_atm[atmname] = HB_atm.pop(retrieve_key(HB_atm,atmcoor))
#            #print(atmname,retrieve_key(HB_atm,atmcoor))
#        if atmname[0] == 'C':
#            if get_str(retrieve_key(CO,atmcoor)) != '':
#                print(atmname,retrieve_key(CO,atmcoor))
#                CO[atmname] = CO.pop(retrieve_key(CO,atmcoor))
##            
##            elif get_str(retrieve_key(Alip_C_BB,atmcoor)) != '':
##                Alip_C_BB[atmname] = Alip_C_BB.pop(retrieve_key(Alip_C_BB,atmcoor))
##            else:
##                Alip_C_Exp[atmname] = Alip_C_Exp.pop(retrieve_key(Alip_C_Exp,atmcoor))
##        if atmname[0] == 'H':
##            NH[atmname] = NH.pop(retrieve_key(NH,atmcoor))
#    
#    print(Alip_C_BB)
#    print(Alip_C_Exp)
#    print(CO)
#    print(NH)
#    CA_list = {'backbone':[],'exposed':[]}
#    for BBkey in Alip_C_BB.keys():
#        CA_list['backbone'].append(BBkey)
#    for Expkey in Alip_C_Exp.keys():
#        CA_list['exposed'].append(Expkey)
#    print(CA_list)
#    
    
def permutation(limitations:list,change_range:dict):#lim_idx:[start,end]
    #input test
    step ={0:1,1:1,2:1,3:0.1}#The HP/HB atm limit must be int, will define expoential coef later
    og_range = {0:[2,8],1:[2,8],2:[2,6],3:[0.2,2]}
    if len(limitations) < len(change_range):
        print('Too much limitations!')
        return 0 
    for keys,each in change_range.items():
        if len(each) > 2:
            print('Wrong format!')
            return 0
        if (each[0] < og_range[keys][0]) or (each[1] > og_range[keys][1]):
            print('Range is too large!')
            return 0
                
    change_idx = list(change_range.keys())
    
    new_range = {}
    for n in change_range.keys(): 
        new_range[n] = []
        k = 0
        while (change_range[n][0]+k*step[n]) <= change_range[n][1]:
            new_range[n].append(round(change_range[n][0]+k*step[n],1))
            k += 1
    print(new_range)
    
    new_lim = list(itertools.product(*new_range.values()))
    print(new_lim)
    print(change_idx)
    all_subset = []
    for each in new_lim:
        temp = list(limitations)
        for n in range(0,len(change_idx)):
            temp[change_idx[n]] = each[n]
        print(temp)
        all_subset.append(temp)
                    
    return all_subset
    

        

    
def get_report(struc,sub_folder,keyword = 'Lig_Only_3_layers.pdb',HP_lim_bb =4,HP_lim_ep = 8, Wat_HB_lim = 3, HP_Crit = 3.5, HB_Crit = 3,HP_Falloff = 0.8, HB_Falloff = 0.6):
    #HB_atm_pair and CA_list stored within a class    
    for file in os.listdir(sub_folder):
        if file.endswith(keyword) and 'Prot' not in file:
            lines = open(sub_folder+'\\'+file,'r').readlines()
            residue, wat = {},{}
            for atm in lines:
                split = atm.split()
                if split[3] == 'WAT' and split[2] == 'O':
                    wat[split[4]+split[5]] = [float(split[6]),float(split[7]),float(split[8])]
                else:
                    try:
                        residue[split[2]] = [float(split[6]),float(split[7]),float(split[8])]
                    except:
                        #print('Error occurred at {}. \n Line as shown: {}'.format(file, atm))
                        pass
            
            #print(residue)
            HB_atm = []
            AliP_C = []
            for key in residue.keys():           
                if key in list(struc.HB_pair.keys()):
                    HB_atm.append(key)
                if key in struc.CA['backbone'] or key in struc.CA['exposed']:
                    AliP_C.append(key)
            #print(HB_atm,AliP_C)   
            
            
            #Hydrophobic
            HP = {}
            for atm in AliP_C:
                HP[atm] = []
                for keys in wat.keys():
                    if HP_HB_Interaction(residue[atm],wat[keys],HP_Crit, HP_Falloff):
                        HP[atm].extend([keys, HP_HB_Interaction(residue[atm],wat[keys],HP_Crit,HP_Falloff)])
                        
                if HP[atm] == []:
                    del HP[atm]
                
            #print('Hydrophobic interactions as shown:',HP)                
            #Intra H bond
            Intra_HB = {}
            O_clash = {}
            for l in range(0,len(HB_atm) -1):
                Intra_HB[HB_atm[l]] = []
                for n in range(l+1, len(HB_atm)):
                    Intra_Crit = Intra_HB_Crit(HB_atm[l],HB_atm[n])
                    if  Intra_Crit == '0':
                        if OO_Clash(residue[HB_atm[l]], residue[HB_atm[n]]):
                            O_clash[HB_atm[l]] = [HB_atm[n],OO_Clash(residue[HB_atm[l]], residue[HB_atm[n]])]
                        continue
                    elif HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit,HB_Falloff):
                        
                        if Intra_Crit == '1':#l is N, n can be both
                            vec_1 = np.array(residue[struc.HB_pair[HB_atm[l]]])-np.array(residue[HB_atm[l]])                        
                            vec_2 = np.array(residue[struc.HB_pair[HB_atm[l]]])-np.array(residue[HB_atm[n]])                        
                            Intra_HB[HB_atm[l]].extend([HB_atm[n],HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit,HB_Falloff),Intra_HB_Geo(vec_1,vec_2)])
                            
#                        #l is O, n is N
                        elif  Intra_Crit == '2':
                            vec1 = np.array(residue[struc.HB_pair[HB_atm[n]]])-np.array(residue[HB_atm[n]])
                            vec2 = np.array(residue[struc.HB_pair[HB_atm[n]]])-np.array(residue[HB_atm[l]])
                            Intra_HB[HB_atm[l]].extend([HB_atm[n],HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit,HB_Falloff),Intra_HB_Geo(vec1,vec2)]) 
                        
                if Intra_HB[HB_atm[l]] == []:
                    del Intra_HB[HB_atm[l]]
            #print(Intra_HB)            
            #Water H bond
            Wat_HB = {}
            
            for atm in HB_atm:
                Wat_HB[atm] = []
                for keys in wat.keys():
                    if HP_HB_Interaction(residue[atm],wat[keys],HB_Crit,HB_Falloff):
                        Wat_HB[atm].extend([keys, HP_HB_Interaction(residue[atm],wat[keys],HB_Crit,HB_Falloff)])
                if Wat_HB[atm] == []:
                    del Wat_HB[atm]
            
            #print(Wat_HB)
            
            #Get score
            score = {'Stability':0,'HP_Ct':0,'HB_Ct':0,'Intra_HBct':0,'O_clash':0}
            penalty_atm ={}
            for keys,values in Wat_HB.items():
                temp = 0
                
                for n in range(0,int(len(values)/2)):
                    temp += (HB_Crit/float(values[2*n +1]))**2
                score['HB_Ct'] += min(temp,Wat_HB_lim)
#                except:
#                    print(keys, values)
#                    return 0
#                score['HB_Ct'] += min((len(values)/2),Wat_HB_lim)
                if len(values) < (2*Wat_HB_lim):
                    penalty_atm[keys] =  Wat_HB_lim -temp
                    #penalty_atm[keys] =  Wat_HB_lim - (len(values)/2)
            #Can be modified for different structure
            
            
            for keys,values in HP.items():
                if keys in struc.CA['backbone']:
                    HP_temp = 0
                    for h in range(0,int(len(values)/2)):
                        HP_temp += (HP_Crit/float(values[2*h +1]))**2
                    score['HP_Ct'] += min(HP_temp,HP_lim_bb)
                    #score['HP_Ct'] += min((len(values)/2),HP_lim_bb)
                if keys in struc.CA['exposed']:   
                    HP_temp2 =0
                    for v in range(0,int(len(values)/2)):
                        HP_temp2 += (HP_Crit/float(values[2*v +1]))**2
                    score['HP_Ct'] += min(HP_temp2,HP_lim_ep)
                    #score['HP_Ct'] += min((len(values)/2),HP_lim_ep)
            for keys,values in Intra_HB.items():
                #We need to debug the intra_HB calculation and add a stability bonus for each GOOD
                #intraHB. 
                score['Intra_HBct'] += (len(values)/3)
                score['Stability'] -= (len(values)/3)
                for atm in list(penalty_atm.keys()):
                    if keys == atm:
                        del penalty_atm[atm]
                    elif atm in values:
                        del penalty_atm[atm]        
            
            for O_values in O_clash.values():
                for l in range(0,int(len(O_values)/2)):
                    score['Stability'] += (3/float(O_values[2*l +1]))**2
                
                score['O_clash'] += (len(O_values)/2)
                
            for each in penalty_atm.values():
                score['Stability'] += each
            
            score['HP percentage'] = score['HP_Ct'] * 100 /(len(struc.CA['exposed'])*HP_lim_ep+len(struc.CA['backbone'])*HP_lim_bb)
            score['Water HB percentage'] = score['HB_Ct']*100 / (len(struc.HB_pair)*Wat_HB_lim)
            
            
    return HP,Intra_HB,Wat_HB, O_clash, score
        
        
def output_report(struc, keyword = 'Pose', output_name = 'Primitive Water Report.txt',opt = 'csv',exp_coef = 0.05):#Do a graph
    os.chdir(struc.Master_folder)
    if opt == 'txt':
        total_ct = open('Count_Stats.txt','w')
    if opt == 'csv':
        total_ct = open('Count_Stats.csv',mode = 'w',newline = '')
        fieldnames = ['Conf', 'HP Interaction', 'Intra HB', 'Wat HB','HP percentage',"Water HB percentage",
                      'hydrophobicity','stability score', 'Estimated population']
        writer = csv.DictWriter(total_ct, fieldnames=fieldnames)
        writer.writeheader()
        
        emptyCt = []
    
    HB_header = list(struc.HB_pair.keys())
    HB_all = {}
    for each in HB_header:
        HB_all[each] ={}
        for atm in HB_header:
            HB_all[each][atm] = 0
            
    for pose in os.listdir(struc.Master_folder):
        if pose.startswith(keyword):
            if os.stat(pose+'\\'+ pose +'_Lig_Only_3_layers.pdb').st_size != 0:            
                try:
                    HP, Intra_HB, Wat_HB, O_clash,score = get_report(struc,pose)
                except Exception as e:
                    print('File {} not processed!\n'.format(pose))
                    print(e)
                    continue

            else:
                emptyCt.append(pose)
                continue

            
            #Count all HB
            for donor,accept in Intra_HB.items():
                try:
                    HB_all[donor][accept[0]] += 1
                except Exception as e:
                    print(e)
                    
            
            
            
            output = open(os.path.join(pose,output_name),'w')
            #Measure of stability

            
            output.writelines('Hydrophobic interactions as shown:\n')
            output.writelines('Alip C\twat_num\tdist\n')           
            for keys,values in HP.items():
                output.write(keys+'\t')
                output.write('\t'.join(map(str,values)))
                output.write('\n')
                
                
            output.writelines('Water-peptide hydrogen bond as shown:\n')
            output.writelines('PepAtm\twat_num\tdist\n')
            
            for keys,values in Wat_HB.items():
                output.write(keys+'\t')
                output.write('\t'.join(map(str,values)))
                output.write('\n')
                
            #print(pose, penalty_atm)    
            output.writelines('Intramolecular hydrogen bonds as shown:\n')
            output.writelines('Don\tAccep\tdist\tAlign\n')
            for keys,values in Intra_HB.items():
                output.write(keys+'\t')
                output.write('\t'.join(map(str,values)))
                output.write('\n')
                
                
                
            output.writelines('Intramolecular Oxygen clash as shown:\n')
            output.writelines('O1\tO2\tdist\n')
            for keys,values in O_clash.items():
                output.write(keys+'\t')
                output.write('\t'.join(map(str,values)))
                output.write('\n')
        
                

            #print('Current score:',score['Stability'])
            output.writelines("Total number of HP interactions: {}\n".format(score['HP_Ct']))
            output.writelines("Total number of intramolecular HB interactions: {}\n".format(score['Intra_HBct']))
            output.writelines("Total number of peptide-water HB interactions: {}\n".format(score['HB_Ct']))
            output.writelines("Total number of intramolecular Oxygen clash: {}\n".format(score['O_clash']))
            output.writelines("Instability score of conformation: {}\n".format(score['Stability']))
            output.close()
            #
#            contribution1 = round(Wat_HBct/0.44,2)
#            contribution2 = contribution1+round(Intra_HBct/0.44,2)
#            contribution3 = round(HPct/0.78,2)
            essential = {'Conf':pose,'HP Interaction':score['HP_Ct'], 'Intra HB': score['Intra_HBct'], 'Wat HB':score['HB_Ct'], 'HP percentage':score['HP percentage'] ,
                         "Water HB percentage":score['Water HB percentage'], 'hydrophobicity': score['HP percentage'] - score['Water HB percentage'],#HP-HB, always
                         'stability score':score['Stability'], 'Estimated population': round(math.exp(-score['Stability']*exp_coef),4)}
            if opt =='txt':
                total_ct.writelines('{}\n'.format(pose))
                total_ct.writelines("HP interactions: {}\n".format(score['HP_Ct']))
                total_ct.writelines("intramolecular HB interactions: {}\n".format(score['Intra_HBct']))
                total_ct.writelines("peptide-water HB interactions: {}\n".format(score['HB_Ct']))
                total_ct.writelines("Total number of intramolecular Oxygen clash: {}\n".format(score['O_clash']))
                total_ct.writelines("Water HB contribution percentage {}%% \n".format(score['Water HB percentage']))
                

            if opt == 'csv':
                writer.writerow(essential)
    total_ct.close()
    
    if len(emptyCt) > 0:
        print('The following {} conformations cannot be processed by Watgen solvation:\n'.format(len(emptyCt)))
        emptyFi = open('emptyFnams.txt','w')
        for emp in emptyCt:
            emptyFi.writelines(emp+'\n')
        emptyFi.close()
    
    #Get a heatmap(triangular) now
    dataset = []
    for line in HB_all.values():
        temp = []
        for instance in line.values():
            temp.append(instance)
        dataset.append(temp)    
    mask =  np.tri(len(HB_header), k=-1)
    print(mask)
    A = np.ma.array(dataset, mask=mask) 
    print(A)
    fig = PLT.figure()
    ax1 = fig.add_subplot(111)
    cmap = CM.get_cmap('coolwarm', 10) # jet doesn't have white color
    cmap.set_bad('w') # default value is 'k'
    img = ax1.imshow(A, interpolation="nearest", cmap=cmap)
    ax1.set_xticks(np.arange(len(HB_header)))
    ax1.set_yticks(np.arange(len(HB_header)))
    ax1.set_xticklabels(HB_header)
    ax1.set_yticklabels(HB_header)
    b, t = PLT.ylim() # discover the values for bottom and top
    b += 0.5 # Add 0.5 to the bottom
    t -= 0.5 # Subtract 0.5 from the top
    PLT.ylim(b, t) # update the ylim(bottom, top) values
     # ta-da!
    ax1.grid(True)

    PLT.colorbar(img,label="HB Count", orientation="vertical",shrink = 0.5)
    PLT.rcParams["figure.figsize"] = (7,7)
    PLT.show()
    new_heatmap = sns.heatmap(A, mask=mask, cmap=cmap, center=0,xticklabels=HB_header, yticklabels=HB_header,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})
    new_heatmap.set(ylim=(b, t))
    
def graph_output(struc,lims,keyword = 'Pose'): #graph multiple results with different interaction lims and exp factorials, so we can predict logP with those parameters. 
    os.chdir(struc.Master_folder)
    
    prediction = {}
    HP_Tol = {}
    P_coef = {}
    
    #Record calculated data
    record = open('LogP fitting records.txt','a')
    record.writelines(starttime+'\n')
    
    score_distri = {}
    for n in range(0,len(lims)):
        ea_distri = {}
        for k in range(0,11):
            ea_distri[k] = 0
        
        empty_Ct = 0
        each_lim = lims[n]
        HP_Tol[n] = {}
        P_coef[n] = 0
        HP_lim_bb,HP_lim_ep,HB_lim, exp_fac = each_lim[0],each_lim[1],each_lim[2],each_lim[3]
        for pose in os.listdir(struc.Master_folder):                
            if pose.startswith(keyword):
                if os.stat(pose+'\\'+ pose +'_Lig_Only_3_layers.pdb').st_size != 0:            
                    try:
                        HP, Intra_HB, Wat_HB, O_clash,score = get_report(struc,pose,HP_lim_bb =HP_lim_bb,HP_lim_ep = HP_lim_ep, Wat_HB_lim = HB_lim)
                        
                        population = max(math.exp(-score['Stability']*exp_fac),10**-7)
                        
                        P_coef[n] += population
                            
                        HP_para = score['HP percentage'] - score['Water HB percentage']
                        #print(HP_para)
                        
                        if HP_para > 0:
                            Log_HP = math.log10(HP_para)
                        elif HP_para < 0:
                            if HP_para < -1:
                                Log_HP = math.log10(abs(HP_para))*-1
                            else:
                                Log_HP = math.log10(abs(HP_para))
                        else:
                            Log_HP = 0
                        Log_HP = Log_HP*3 -1     
                        HP_Tol[n][pose] = [population, Log_HP]
                        if score['Stability'] <= 10:
                            ea_distri[score['Stability']] += 1   
                        
                    except:
                        #print('File {} not processed!\n'.format(pose))
                        
                        continue
                else:
                    #print('{} is empty'.format(pose))
                    empty_Ct +=1
                    continue
        score_distri[n] = ea_distri        
        
            
        print('{} files are found to be empty in this run.'.format(empty_Ct))
        #Get population(exp factorial)+add up population            
        #Sum of log HP--treat negative values, this is predicted 
        #Add results and plot
#    print(HP_Tol)
    print(P_coef)
    for n in HP_Tol.keys():
        logP_Tol = 0
        for conf in HP_Tol[n].values():
            each = conf[0] *conf[1]/P_coef[n]
            logP_Tol+= each
        
        prediction[n] = logP_Tol
        record.writelines('Limitations as shown: backbone HP limit {}, HP limit {}, HB limit {}, exponential factor {} \n'.format(lims[n][0],lims[n][1],lims[n][2],lims[n][3])) 
        
        for keys,items in score_distri[n].items():
            record.writelines('{} conformers scored {}.\n'.format(items,keys))
        
        
        record.writelines('Prediction value of logP is: {}\n\n'.format(logP_Tol))
        
        
        #Save plot
    
    print(prediction)
#    x,y, a,b = [],[],[],[]
#    for n in range(0,3):
#        x.append(lims[n][0])
#        a.append(prediction[n])
#    for n in range(3,len(lims)):
#        y.append(lims[n][0])
#        b.append(prediction[n])
#    
#        
#    plt.plot(x, a, 'o',)
#    plt.plot(y, b, 'o',)
#    plt.show()
#        

#Testing block here
#if __name__ == '__main__':
#    pre_add = 'C:\\People\Kaichen\\PEPEDIT\\A6-2NMe_2DAA\\poses_m\\out000001m.pdb'
#    solv_add = 'C:\\Users\\Haworth_Lab2\\Desktop\\Solvate\\Solvate\\Inputs\\A6-2NMe_2DAA\\Pose_1\\Pose_1_compare_waters.pdb'
#    trace_rename_atm(pre_add,solv_add)    
