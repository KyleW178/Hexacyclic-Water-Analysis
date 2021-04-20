# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 23:32:55 2021

@author: Haworth_Lab2
"""

import os,math,csv,pytz,itertools
import numpy as np
import matplotlib.pyplot as plt
import datetime
        
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

def OO_Clash(c1,c2,crit = 4):
    distance = math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)
    if distance < crit:
        return str(round(distance, 3))
    else:
        return False

#Why this falloff design doesnot work???
#def HP_HB_Interaction(coo1,coo2,crit):
#    distance = math.sqrt((coo1[0]-coo2[0])**2 + (coo1[1]-coo2[1])**2 + (coo1[2]-coo2[2])**2)
#    falloff = 0.7
#    if (crit - falloff)<= distance <= (crit + falloff):#HP crit 3.5, HB crit 3
#        return str(round(distance, 3))
#    else:
#        return False

def HP_HB_Interaction(coo1,coo2,crit):
    distance = math.sqrt((coo1[0]-coo2[0])**2 + (coo1[1]-coo2[1])**2 + (coo1[2]-coo2[2])**2)
    if distance <= crit:#HP crit 3.5, HB crit 3
        return str(round(distance, 3))
    else:
        return False

#Name: N01, O04, so on
#Please make a new one... this does not work
#Consider: regular+ O-O clash(return 0)
def Intra_HB_Crit(name1, name2, aa = 6):
    a, b = get_num(name1), get_num(name2)
    if abs(a-b) in range(2, aa-1):
        if (name1[0]  == 'O' and name2[0] == 'O'):
            return '0'
        else:
            if name1[0] == 'N':
                return 'x'
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
    #print('Angle is:',ang)
    if  final > limit:
        return 'NG'
    else:
        return final
    

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
    

        

    
def get_report(struc,sub_folder,keyword = 'Lig_Only_3_layers.pdb',HP_lim_bb =4,HP_lim_ep = 7, Wat_HB_lim = 3, HP_Crit = 3.5, HB_Crit = 3):
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
#                        print('Error occurred at {}. \n Line as shown: {}'.format(file, atm))
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
                    if HP_HB_Interaction(residue[atm],wat[keys],HP_Crit):
                        HP[atm].extend([keys, HP_HB_Interaction(residue[atm],wat[keys],HP_Crit)])
                if HP[atm] == []:
                    del HP[atm]
                
            #print('Hydrophobic interactions as shown:',HP)                
            #Intra H bond
            Intra_HB = {}
            O_clash = {}
            for l in range(0,len(HB_atm) -1):
                Intra_HB[HB_atm[l]] = []
                for n in range(l+1, len(HB_atm)):
                    
                    if HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit):
                        Intra_Crit = Intra_HB_Crit(HB_atm[l],HB_atm[n])
                        if  Intra_Crit == 'x':#l is N, n can be both
                            try:
                                vec_1 = np.array(residue[struc.HB_pair[HB_atm[l]]])-np.array(residue[HB_atm[l]])
                            
                                vec_2 = np.array(residue[struc.HB_pair[HB_atm[l]]])-np.array(residue[HB_atm[n]])
                            
                                Intra_HB[HB_atm[l]].extend([HB_atm[n],HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit),Intra_HB_Geo(vec_1,vec_2)])
                            except:
                                print('Problem occurred with fetching intra HB!\n')
#                        #l is O, n is N
                        elif  Intra_Crit == '2':
                            vec1 = np.array(residue[struc.HB_pair[HB_atm[n]]])-np.array(residue[HB_atm[n]])
                            vec2 = np.array(residue[struc.HB_pair[HB_atm[n]]])-np.array(residue[HB_atm[l]])
                            Intra_HB[HB_atm[l]].extend([HB_atm[n],HP_HB_Interaction(residue[HB_atm[l]], residue[HB_atm[n]], HB_Crit),Intra_HB_Geo(vec1,vec2)]) 
                        elif  Intra_Crit == '0':
                            O_clash[HB_atm[l]] = [HB_atm[n],OO_Clash(residue[HB_atm[l]], residue[HB_atm[n]])]       
    #                            except:
#                                print('Error at {} and {}'.format(HB_atm[l],HB_atm[n]))
                if Intra_HB[HB_atm[l]] == []:
                    del Intra_HB[HB_atm[l]]
            #print(Intra_HB)            
            #Water H bond
            Wat_HB = {}
            
            for atm in HB_atm:
                Wat_HB[atm] = []
                for keys in wat.keys():
                    if HP_HB_Interaction(residue[atm],wat[keys],HB_Crit):
                        Wat_HB[atm].extend([keys, HP_HB_Interaction(residue[atm],wat[keys],HB_Crit)])
                if Wat_HB[atm] == []:
                    del Wat_HB[atm]
            
            #print(Wat_HB)
            
            #Get score
            score = {'Stability':0,'HP_Ct':0,'HB_Ct':0,'Intra_HBct':0,'O_clash':0}
            penalty_atm ={}
            for keys,values in Wat_HB.items():
                score['HB_Ct'] += min((len(values)/2),Wat_HB_lim)
                if len(values) < (2*Wat_HB_lim):
                    penalty_atm[keys] =  Wat_HB_lim - (len(values)/2)
            #Can be modified for different structure
            for keys,values in HP.items():
                if keys in struc.CA['backbone']:
                    score['HP_Ct'] += min((len(values)/2),HP_lim_bb)
                if keys in struc.CA['exposed']:   
                    score['HP_Ct'] += min((len(values)/2),HP_lim_ep)
            for keys,values in Intra_HB.items():
                score['Intra_HBct'] += (len(values)/3)
                for atm in list(penalty_atm.keys()):
                    if keys == atm:
                        del penalty_atm[atm]
                    elif atm in values:
                        del penalty_atm[atm]        
            
            for values in O_clash.items():
                
                score['Stability'] += (len(values)/2)*3
                
                score['O_clash'] += (len(values)/2)
                
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
    for pose in os.listdir(struc.Master_folder):
        if pose.startswith(keyword):
            if os.stat(pose+'\\'+ pose +'_Lig_Only_3_layers.pdb').st_size != 0:            
                try:
                    HP, Intra_HB, Wat_HB, O_clash,score = get_report(struc,pose)
                except Exception as e:
                    print('File {} not processed!\n'.format(pose))
                    print(e)
                    #continue
                    return 0 
            else:
                emptyCt.append(pose)
                continue
#            HPct, Intra_HBct, Wat_HBct, O_clashct =0,0,0,0
            
            
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
        print(n)
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

    