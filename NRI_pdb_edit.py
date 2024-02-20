# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 16:40:41 2023

@author: Ehsan
"""

import pandas as pd
import os
import glob
import re

# path= '/media/arma/DATA/S-Ehsan/Automate_sMD/HIV_cocrystal_inhibitors/prepared_molecules/'
path= 'C:/Users/ehsan/Desktop/'
os.chdir(path)
files= glob.glob('*.pdb')

def NRI_pdb (file):
    with open ('{}.pdb'.format(file),'r') as pdb:
        pdbtxt=pdb.read()
    pdbedit=pdbtxt.replace(' X ', '  ')  
    pdbedit=pdbedit.replace(' B ', '  ')
    pdbedit=pdbedit.replace(' A ', '  ') 
    pdbedit=pdbedit.replace('   C2', '     ') 
    with open('{}-first.pdb'.format(file), 'w') as pdb_file:
        pdb_file.write(pdbedit)
    with open ('{}-first.pdb'.format(file),'r') as pdb:
        line = pdb.readlines()
    df = pd.DataFrame(line)
    mask = (df[0].str.startswith('ATOM')) & (~df[0].str.contains('UNK')) & (~df[0].str.contains('CA'))
    # Apply the mask to filter the DataFrame
    filtered_df = df[~mask]
    ddf = filtered_df.reset_index(drop=True)
    contains_end = ddf.apply(lambda row: any('END' in cell for cell in row), axis=1)
    contains_pseu = ddf.apply(lambda row: any('pseu' in cell for cell in row), axis=1)
    need_index = []
    for i in range(len(contains_end)):
        if contains_end[i] == True:
            ddf.loc[i,:] = 'ENDMDL'
            need_index.append(i)
        if contains_pseu[i] == True:
            ddf.drop(i, inplace = True)
    rows=[]
    index = []
    for z in range (2,2003):
        new_rows= {line[0].replace('\n',''):'MODEL     {}'.format(z)}
        rows.append(new_rows)
    for x in range(len(rows)):
        ddf.loc[need_index[x]+0.5,:]= rows[x].get(line[0].replace('\n',''))
    ddf = ddf.sort_index().reset_index(drop=True)
    dff = pd.DataFrame(ddf)
    dff.to_csv('{}-sec.pdb'.format(file),header=None ,index=None)
    with open ('{}-sec.pdb'.format(file),'r') as pdb:
        lines = [line.strip('"\n') for line in pdb]
    line = [line.strip() for line in lines if line.strip() != '']
    line[0]='MODEL     1'
    renumbered_data = []
    model_number = 0
    atom_count = 0
    for l in line:
        if l.startswith('MODEL'):
            model_number = int(l.split()[1])
            atom_count = 0
            renumbered_data.append(l)
        elif l.startswith('ATOM'):
            atom_count += 1
            line_parts = l.split()
            line_parts[1] = str(atom_count)
            renumbered_line = ' '.join(line_parts)
            renumbered_data.append(renumbered_line)
        elif l.startswith('ENDMDL'):
            renumbered_data.append(l)
    for i in range(len(renumbered_data)):
        if len(renumbered_data[i]) > 40:
            if len(renumbered_data[i].split()[1]) == 1 :
                split_data = renumbered_data[i].split()
                split_data[4] = split_data[1]
                renumbered_data[i] = ' '.join(split_data)
                # renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[4], renumbered_data[i].split()[1])
            if len(renumbered_data[i].split()[1]) == 2 :
                split_data = renumbered_data[i].split()
                split_data[4] = split_data[1]
                renumbered_data[i] = ' '.join(split_data)
                #renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[4], renumbered_data[i].split()[1])
            if len(renumbered_data[i].split()[1]) == 3 :
                # if int(renumbered_data[i].split()[1]) > 300:
                #     if int(renumbered_data[i].split()[1]) == 332:
                if renumbered_data[i].split()[3] == 'UNK':
                    split_data = renumbered_data[i].split()
                    split_data[4] = split_data[1]
                    renumbered_data[i] = ' '.join(split_data)
#                    renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[4], renumbered_data[i].split()[1])
                if renumbered_data[i].split()[2] == 'UNK':
                    split_data = renumbered_data[i].split()
                    split_data[3] = split_data[1]
                    renumbered_data[i] = ' '.join(split_data)
                    #renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[3], renumbered_data[i].split()[1])
                else:
                    split_data = renumbered_data[i].split()
                    split_data[4] = split_data[1]
                    renumbered_data[i] = ' '.join(split_data)
                    #renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[4], renumbered_data[i].split()[1])
            else:
                split_data = renumbered_data[i].split()
                split_data[4] = split_data[1]
                renumbered_data[i] = ' '.join(split_data)
                #renumbered_data[i] = renumbered_data[i].replace(renumbered_data[i].split()[4], renumbered_data[i].split()[1])
    line = renumbered_data
    new_line = '\n'.join(line)   
    lines = new_line.split('\n')
    output_file_path ='{}-1.pdb'.format(file)
    # Open the output file in write mode and write the modified lines
    with open(output_file_path, 'w') as output_file:
        output_file.write(new_line)
    return print('Pdb File of {} For NRI is Ready!!'.format(file))

for n in files:
    file=n.split('.')[0]
    NRI_pdb(file)
