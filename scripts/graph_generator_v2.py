#!/usr/bin/python

import os
import os.path
import argparse
#import re
#from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
import pygraphviz as pgv
import math
from solution_parser import extract_from_file, sspace_extract


datatypes = ['inpt', 'expt', 'whpm', 'dist', 'heur', 'brpu', 'flow', 'sspa']

#GRAPH CREATION METHOD
def create_graph_1o1n(datafile, datatype, solarray):
    dtype_2_glabel = {'inpt':"Input data of ", 'expt':"Expected solution of "
                    , 'whpm':"Weighted path model solution of "
                    , 'dist':"Distance based model solution of "
                    , 'heur':"Genetic algorithm model solution of "
                    , 'brpu':"Branch and prune model solution of "
                    , 'flow':"Flow model solution of "
                    , 'sspa':"Sppace solution of "}
    dtype_2_gname = {'inpt':'INPT', 'expt':"EXPT", 'whpm':"WHPM", 'dist':"DIST"
                   ,'heur':"HEUR", 'brpu':"BRPU", 'flow':"FLOW", 'sspa':'SSPACE'}
    cord_2_color = {'R':'indianred1', 'F':'mediumslateblue', 'O':'limegreen', 'I':'gray'}
    o1o2_2_color = {'RR':'indianred1', 'FF':'mediumslateblue'
                  ,'RF':'indianred1;0.5:mediumslateblue'
                  , 'FR':'mediumslateblue;0.5:indianred1'
                  , 'OF':'limegreen;0.5:mediumslateblue' 
                  , 'OR':'limegreen;0.5:indianred1'
                  , 'FO':'mediumslateblue;0.5:limegreen'
                  , 'RO':'indianred1;0.5:limegreen'}

    nbsol = len(solarray)
    dataname = os.path.basename(datafile)
    foldername = datafile.split('.')[0]
    try:
        foldername = foldername.split('_',1)[1]
    except IndexError:
        foldername = foldername

    for solidx in range (0, nbsol):
        list_unitigs, dict_unitg2len, dict_unitg2cov, list_links =  solarray[solidx][:]

        if datatype != 'inpt':
            G=pgv.AGraph(strict = False, directed = True)
            for unitig in list_unitigs:
                unr = unitig
                try:
                    ulen = dict_unitg2len[unr.split('-')[0]]
                except KeyError:
                    ulen = '?'
                try:
                    ucov = dict_unitg2cov[unr.split('-')[0]]
                except KeyError:
                    ucov = '?'

                if datatype == 'whpm':
                    G.add_node(unr, label="{ unitig " + unr.split('-')[0] +"| occ. "+ unr.split('-')[1]+" of "+ str(ucov) +"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')
                if datatype == 'dist':
                    G.add_node(unr, label="{ unitig " + unr.split('-')[0] +"| occ. "+ unr.split('-')[1]+" of "+ str(ucov) +"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')
                if datatype == 'flow':
                    G.add_node(unr, label="{ unitig " + unr.split('-')[0] +"| len. "+ str(ulen) +"| occ. "+ unr.split('-')[1]+" of "+ str(ucov) +"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')
                if datatype == 'expt':
                    G.add_node(unr, label="{ unitig " + unr.split('-')[0] +"| len. "+ str(ulen) +"| occ. "+ unr.split('-')[1]+" of "+ str(ucov) +"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')
                if datatype == 'sspa':
                    G.add_node(unr, label="{ unitig " + unr+"| len. "+ str(ulen)+"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')                               

            for link in list_links:
                unr1, unr2, u_ori1, u_ori2, uudist = link[:]
                G.add_edge(unr1, unr2, label=uudist, style='bold', penwidth=3, color=o1o2_2_color[u_ori1+u_ori2], fontname='AvantGarde-Book', fontsize="10.0")
                nunr1=G.get_node(unr1)
                nunr2=G.get_node(unr2)
                nunr1.attr['color'] = cord_2_color[u_ori1]
                nunr2.attr['color'] = cord_2_color[u_ori2]


            G.graph_attr['fontname'] = 'AvantGarde-Book'
            G.graph_attr['label'] = dtype_2_glabel[datatype] + dataname + " (1 occurrence, 1 node graph), Solution #" + str(solidx)
            G.layout('circo')
            #G.write(foldername + "/" + datafile + '_G_' + dtype_2_gname[datatype] + ".dot")
            G.draw(foldername + "/" + dataname + '_G_' + dtype_2_gname[datatype] +"_"+str(solidx) + "_.png")
        
        elif datatype == 'inpt':
            G=pgv.AGraph(strict = False, directed = True) # graphical representation of literal content of input file
            G2=pgv.AGraph(strict = False, directed = True) # graphical representation of the modeled problem

            for unitig in list_unitigs:
                unr = unitig
                ulen = dict_unitg2len[unr.split('-')[0]]
                ucov = dict_unitg2cov[unr.split('-')[0]]

                G.add_node(unr, label="{ unitig " + unr +"| cov. "+ str(ucov) +"| len. "+ str(ulen) +"}", shape = 'Mrecord', \
                               style='filled', fontname='AvantGarde-Book', fontcolor='black')

                ucov_max = int(ucov[1])
                for uidx in range (1, ucov_max+1):
                    unruocc_F = unr+"-"+str(uidx)+"-F"
                    unruocc_R = unr+"-"+str(uidx)+"-R"
                    G2.add_node(unruocc_F, label="{ unitig " + unr +"| occ. "+ str(uidx) +" of "+str(ucov_max) +"| len. "+ str(ulen) +"}", shape = 'Mrecord', \
                                   style='filled',  color=cord_2_color['F'], fontname='AvantGarde-Book', fontcolor='black')                    
                    G2.add_node(unruocc_R, label="{ unitig " + unr +"| occ. "+ str(uidx) +" of "+str(ucov_max) +"| len. "+ str(ulen) +"}", shape = 'Mrecord', \
                                   style='filled',  color=cord_2_color['R'], fontname='AvantGarde-Book', fontcolor='black')              
            for link in list_links:
                unr1, unr2, u_ori1, u_ori2, uudist = link[:]
                G.add_edge(unr1, unr2, label=uudist, style='bold', penwidth=3, color=o1o2_2_color[u_ori1+u_ori2], fontname='AvantGarde-Book', fontsize="10.0")
                #nunr1=G.get_node(unr1)
                #nunr2=G.get_node(unr2)
                #nunr1.attr['color'] = cord_2_color[u_ori1]
                #nunr2.attr['color'] = cord_2_color[u_ori2]

                for uidx_1 in range(1, int(dict_unitg2cov[unr1][1])+1):
                    for uidx_2 in range(1, int(dict_unitg2cov[unr2][1])+1):
                        unruocc_1_u_ori1 = unr1+"-"+str(uidx_1)+"-"+u_ori1
                        unruocc_2_u_ori2 = unr2+"-"+str(uidx_2)+"-"+u_ori2
                        G2.add_edge(unruocc_1_u_ori1, unruocc_2_u_ori2, label=uudist, style='bold', penwidth=3, color=o1o2_2_color[u_ori1+u_ori2], fontname='AvantGarde-Book', fontsize="10.0")

            G.graph_attr['fontname'] = 'AvantGarde-Book'
            G.graph_attr['label'] = dtype_2_glabel[datatype] + dataname + " (1 unitig, 1 node graph), Solution #" + str(solidx)
            G.layout('circo')
            #G.write(foldername + "/" + datafile + '_G_' + dtype_2_gname[datatype] + ".dot")
            G.draw(foldername + "/" + dataname + '_G_' + dtype_2_gname[datatype] + "_.png")

            G2.graph_attr['fontname'] = 'AvantGarde-Book'
            G2.graph_attr['label'] = dtype_2_glabel[datatype] + dataname + " (1 unitig occurrence, 1 node graph), Solution #" + str(solidx)
            G2.layout('circo')
            #G.write(foldername + "/" + datafile + '_G_' + dtype_2_gname[datatype] + ".dot")
            G2.draw(foldername + "/" + dataname + '_G2_' + dtype_2_gname[datatype] + "_.png")


    return G

#EXECITIONER METHOD
def executioner(graph, filegraph, evidence):
    dtype_2_msgs = {'inpt':'Input data', 'expt':'Expected solution', \
                  'whpm':'Weighted path model solutions', \
                  'dist':'Distance based model solution', \
                  'brpu':'Branch and prune model solutions', \
                  'heur':'Genetic algorithm model solutions', \
                  'flow':'Flow model solutions', \
                  'sspa':'Sppace solutions'}
    if evidence == 0:
        array_of_multiple_solutions = extract_from_file(graph, filegraph)
        if len(array_of_multiple_solutions) != 0:
            print("Creating graph for the "+dtype_2_msgs[graph]+".")
            create_graph_1o1n(filegraph.split('/')[-1], graph, array_of_multiple_solutions)
        else:
            print("The file contains no solution(s)", filegraph)
    else:
        array_of_multiple_solutions = sspace_extract(filegraph, evidence)
        if len(array_of_multiple_solutions) != 0:
            print("Creating graph for the "+dtype_2_msgs[graph]+".")
            create_graph_1o1n(filegraph.split('/')[-1], graph, array_of_multiple_solutions)
        else:
            print("The file contains no solution(s)", filegraph)        

    return None

#MAIN METHOD
def main():
    parser=argparse.ArgumentParser(description="Creating (genome) graphs")
    parser.add_argument('-g', '--graph', type=str, choices=['dist', 'whpm', 'brpu', 'heur', 'expt', 'inpt', 'flow', 'sspa']
                        , help='choose the types of graphs to generate for all or a single file.')
    parser.add_argument('-f', '--file', type=str, help='path of the file for which the solution graph will be generated')
    parser.add_argument('-e', '--evidence', type=str, help='evidence file for sspace solution')
    args=parser.parse_args()

    pathfile = args.file 
    foldername = (os.path.basename(pathfile)).split('.')[0]
    try:
        foldername = foldername.split('_',1)[1]
    except IndexError:
        foldername = foldername
    
    if os.path.exists(foldername) == False:
        os.mkdir(foldername)
    if args.graph in datatypes:
        if not args.evidence:
            executioner(args.graph,args.file, 0)
        elif args.evidence:
            executioner(args.graph,args.file,args.evidence)
    return args

if __name__ == "__main__":
    args=main()
    print args