#!/usr/bin/python

import argparse
import re
import sys

dict_name2abbr = {
                 'SSPACE':'sspc', 
                 'Input':'inpt', 
                 'Expected':'expt', 
                 'Distance_based':'dist', 
                 'Weighted_path':'whpm', 
                 'Flow_model':'flow',
                 'Sspace scaffolder':'sspa'
                 }


def extract_from_file(gtype, pathfile):
#Comment
    array_of_multiple_solutions = []

    list_unitigs = []
    list_multiplied_unitigs = []
    dict_unitg2len = {}
    dict_unitg2cov = {}

    list_links = []
    

    if gtype == 'inpt':
        re_links = re.compile(r"[0-9]*__len__[0-9]*__[FR] ")
        re_unitg = re.compile(r"[0-9]*__len__[0-9]* [0-9]* [0-9]*")
        f = open(pathfile, 'r')
        for line in f:
            if re_unitg.match(line):
                uname, ulen, ucov_min, ucov_max = line.split()
                uname = uname.split('_')[0]
                list_unitigs.append(uname)
                dict_unitg2len[uname] = ulen
                dict_unitg2cov[uname] = [ucov_min, ucov_max]
            if re_links.match(line):
                if line.split() == 4:
                    u1, u2, d1, d2 = line.split()
                    uname1, _, ulen1, u_ori1 = u1.split('__')
                    uname2, _, ulen2, u_ori2 = u2.split('__')
                    list_links.append([uname1, uname2, u_ori1, u_ori2, str((int(d1)+int(d2))/2)])
                else:
                    u1, u2, dist = line.split()
                    uname1, _, ulen1, u_ori1 = u1.split('__')
                    uname2, _, ulen2, u_ori2 = u2.split('__')
                    list_links.append([uname1, uname2, u_ori1, u_ori2, dist])
        f.close()
        array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])
    
    if gtype == 'expt':
        regex_solline=re.compile(r"[0-9]+\s+[0-9]+\s+[0-9]+__len__[0-9]+\s+[0-9]+\s+[0-9]+\s+[0-9]+\s+(\+|-)")
        f=open(pathfile, 'r')

        pre_link_list = []

        uidx=0
        for line in f:
            if regex_solline.match(line):
                try:
                    cord_g1, cord_g2, unitg, umap, cord_u1, cord_u2, u_ori, uudist = line.split()
                    unr, _, ulen = unitg.split('__')
                    uudist = str((-int(uudist)))
                except ValueError:
                    cord_g1, cord_g2, unitg, umap, cord_u1, cord_u2, u_ori = line.split() ; uudist = 'e'
                    unr, _, ulen = unitg.split('__')
            
                fulen=float(ulen)
                fumap=float(umap)
                if float(fumap/fulen)>0.7:
                    dict_unitg2len[unr] = ulen
                    if unr not in dict_unitg2cov.keys():
                        dict_unitg2cov[unr] = 1
                    else:
                        dict_unitg2cov[unr] = dict_unitg2cov[unr] + 1

                    list_unitigs.append(unr+'-'+str(dict_unitg2cov[unr]))


                    pre_link_list.append(((unr+'-'+str(dict_unitg2cov[unr]) +'_'+ u_ori+'_'+ ulen +'_'+ \
                              umap +'_'+ cord_u1 +'_'+ cord_u2 +'_'+ cord_g1 +'_'+ cord_g2 + '_' + uudist +'_' + str(uidx)) \
                             .replace("_-_", "_R_")).replace("_+_","_F_"))
                    uidx=uidx+1

        list_links=create_links_from_prelinks(pre_link_list)
        f.close()

        array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])
        req_list_links = [reverse_equivalent_link(l) for l in list_links]
        array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, req_list_links])

    if gtype == 'whpm':
        regex_cntg = re.compile(r"[ma]_[0-9]*_[0-9]*_[0-9]*_[FR]=1.0")
        
        f = open(pathfile, 'r')
        nbu = f.readline()

        pre_link_list= [0]*(int(nbu) + 1)

        for line in f:
            if regex_cntg.search(line):
                unitigs = line.split()

                for u in unitigs[:-1]:
                    if '=1.0' in u and "(m_" in u:
                        _, unr, uocc, urank, u_ori = u.split('_')
                        u_ori = u_ori.replace('=1.0)','')
                        uocc = str(int(uocc)+1)

                        if unr not in dict_unitg2cov.keys():
                            dict_unitg2cov[unr] = 1
                        else:
                            dict_unitg2cov[unr] = dict_unitg2cov[unr] + 1
                    
                        list_unitigs.append(unr+'-'+uocc)

                        pre_link_list[int(urank)] = unr+'-'+uocc+'_'+u_ori+'_'+uocc+'_e_e_e_e_e_0_'+urank

            if "END" in line:
                list_unitigs = [x for x in list_unitigs if x !=0]
                pre_link_list = [x for x in pre_link_list if x !=0]
                list_links = create_links_from_prelinks(pre_link_list)
                array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])
                list_unitigs = []; dict_unitg2len = {}; dict_unitg2cov={}; list_links=[]

        f.close()

    if gtype == 'dist' or gtype == 'brpu' or gtype == 'heur':
        regex_dist = re.compile(r"m[0-9]*[F-R][0-9]*")

        unitigs = []
        f = open(pathfile, 'r') 
        for line in f:
            if regex_dist.match(line) or 'ALEXANDRINA END' in line:
                if 'END' not in line:
                    unitigs.append(line)
                elif 'END' in line:
                    uidx=0
                    pre_link_list = []
                    for u in unitigs:
                        k, h = u.split()
                        unr = k[1:-2] ; u_ori = k[-2:-1] ; uudist = float(h)
                        if unr not in dict_unitg2cov.keys():
                            dict_unitg2cov[unr] = 1
                        else:
                            dict_unitg2cov[unr] = dict_unitg2cov[unr] + 1

                        unr = unr+'-'+str(dict_unitg2cov[unr])
                        unitig = unr + '_' + u_ori + '_e_e_e_e_e_e_' + str(uudist) + '_' + str(uidx)
                        list_unitigs.append(unr)
                        pre_link_list.append(unitig)
                        uidx=uidx+1

                    list_links = create_links_from_prelinks(pre_link_list)
                    array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])
       
        f.close()

    if gtype == 'flow':

        f = open(pathfile, 'r')
        
        re_links = re.compile(r"\(\'[0-9]*")
        re_unitg = re.compile(r"C_[0-9]*_[0-9]*_[F-R]\s+[0-9]+")
        re_solution_separator = re.compile(r"Solution\s+[0-9]+")

        for line in f:
            if re_unitg.match(line) or re_links.match(line):
                if re_unitg.match(line):
                    unitig = line.split()[0]
                    _, unr, uocc, u_ori = unitig.split('_')
                    ulen = line.split()[1]

                    dict_unitg2len[unr] = ulen
                    if unr not in dict_unitg2cov.keys():
                        dict_unitg2cov[unr] = 1
                    else:
                        dict_unitg2cov[unr] = dict_unitg2cov[unr] + 1
                    list_unitigs.append(unr+'-'+uocc)
                elif re_links.match(line):
                    unitig = line.split(',')[0]
                    uudist = line.split(',')[1]

                    step, _, _, unr1, uocc1, u_ori1, unr2, uocc2, u_ori2 = unitig.split('_')

                    uudist = uudist[:-2]
                    u_ori2 = u_ori2[:-1]

                    link = [unr1+'-'+uocc1, unr2+'-'+uocc2, u_ori1, u_ori2, uudist]

                    list_links.append(link)
            
            elif re_solution_separator.match(line):
                list_links = []; list_unitigs = []; dict_unitg2cov = {}; dict_unitg2len = {}; pre_link_list = []

                array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])

    return array_of_multiple_solutions
    
def create_links_from_prelinks(pre_link_list):

    list_links = []
    unr1, u_ori1, _, _, _, _, _, _, dist1, uidx1 = pre_link_list[0].split('_')
    unr2, u_ori2, _, _, _, _, _, _, dist2, uidx2 = pre_link_list[-1].split('_')
    link = [unr2, unr1, u_ori2, u_ori1, dist2]
    list_links.append(link)

    len_pll = len(pre_link_list)
    n = 0

    while n < len_pll-1:
        unr1, u_ori1, _, _, _, _, _, _, dist1, uidx1 = pre_link_list[n].split('_')
        unr2, u_ori2, _, _, _, _, _, _, dist2, uidx2 = pre_link_list[n+1].split('_')
        link = [unr1, unr2, u_ori1, u_ori2, dist2]
        n = n+1
        list_links.append(link)

    return list_links

def reverse_equivalent_link(link):

    dict_ori2eqori = {'FF': 'RR', 'RR': 'FF', 'RF':'RF', 'FR':'FR'}
    unr1, unr2, u_ori1, u_ori2, dist = link[:]
    neworis = dict_ori2eqori[u_ori1+u_ori2]
    req_link = [unr2, unr1, neworis[0], neworis[1], dist]

    return req_link


def sspace_extract(formatted_file, evidence_file):
    #formatted_file is "standard_output.extension_evidence.txt" and evidence_file is "standard_output.final.evidence"

    array_of_multiple_solutions = []

    list_unitigs = []
    list_multiplied_unitigs = []
    dict_unitg2len = {}
    dict_unitg2cov = {}

    list_links = []
    pre_link_list = []
    
    dict_sspacenrc_2_unr = {}


    f = open(formatted_file, 'r')
    for line in f:
        if 'contig' in line:
            unr = line.split(':')[1][:-1]
            sspacenrc = (line.split('|')[0]).split('contig')[1]
            unr, _, ulen = unr.split('__')[:]

            dict_sspacenrc_2_unr[sspacenrc] = unr
            list_unitigs.append(unr)
            #dict_unitg2len[unr] = ulen
            dict_unitg2cov[unr] = '?'
    f.close()

    e = open(evidence_file, 'r')
    for line in e:
        arrline = line.split('|')
        if '_tig' in line:
            sspacenrc = arrline[0].split('_tig')[1]
            orientation = arrline[0].split('_tig')[0]
            size = arrline[1][4:]
            dict_unitg2len[dict_sspacenrc_2_unr[sspacenrc]] = size.rstrip()
            try :
                gap = (arrline[3].split('gaps')[1])[:-1]
                pre_link_list.append(dict_sspacenrc_2_unr[sspacenrc]+"_"+orientation+"_e_e_e_e_e_e_"+gap+"_0")
            except IndexError:
                pre_link_list.append(dict_sspacenrc_2_unr[sspacenrc]+"_"+orientation+"_e_e_e_e_e_e_end_0")
    e.close()

    len_pll = len(pre_link_list)
    n = 0
    while n < len_pll-1:
        unr1, u_ori1, _, _, _, _, _, _, dist1, uidx1 = pre_link_list[n].split('_')
        unr2, u_ori2, _, _, _, _, _, _, dist2, uidx2 = pre_link_list[n+1].split('_')
        if dist1 == 'end':
            n = n+1
            continue
        else:
            link = [unr1, unr2, u_ori1.upper(), u_ori2.upper(), dist2]
            n = n+1
            list_links.append(link)



    array_of_multiple_solutions.append([list_unitigs, dict_unitg2len, dict_unitg2cov, list_links])
    return array_of_multiple_solutions


def main():
    parser=argparse.ArgumentParser(description="Creating link and node lists from scaffolding\
                                   solutions, input data or expected solution")
    parser.add_argument('-t', '--type', type=str,
                        choices=dict_name2abbr.values(),
                        help='type of the .txt file to parse')
    parser.add_argument('-f', '--file', type=str, 
                        help='path to .txt file containing the scaffolding solution')
    parser.add_argument('-e', '--evidence', type=str, 
                        help='evidence file for sspace solution')
    args=parser.parse_args()

    if args.type != 'sspa':
        array_of_multiple_solutions = extract_from_file(args.type, args.file)
    elif args.type == 'sspa':
        array_of_multiple_solutions = sspace_extract(args.file, args.evidence)
    print(array_of_multiple_solutions)

if __name__ == "__main__":
  args=main()



    