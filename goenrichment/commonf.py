from rpy2.robjects.vectors import DataFrame, FloatVector, IntVector, StrVector, ListVector
import numpy
from collections import OrderedDict
from goenrichment import gosetting
import rpy2.robjects as robjects
import os
import csv

def recurList(data):
    rDictTypes = [DataFrame, ListVector]
    rArrayTypes = [FloatVector, IntVector]
    rListTypes = [StrVector]
    if type(data) in rDictTypes:
        return OrderedDict(zip(data.names, [recurList(elt) for elt in data]))
    elif type(data) in rListTypes:
        return [recurList(elt) for elt in data]
    elif type(data) in rArrayTypes:
        return numpy.array(data)
    else:
        if hasattr(data, "rclass"):
            raise KeyError('Could not proceed, type {} is not defined'.format(type(data)))
        else:
            return data


def GOEnri(parameters):
    robjects.r('''library(GOstats)''')
    robjects.r('''library(GSEABase)''')
    build_GO = '''BostGO<-read.csv("%s", sep="%s")''' % (os.path.join(gosetting.UPLOAD_FOLDER, parameters['association']), parameters['dassociation'])
    robjects.r(build_GO)
    universe_list = '''universe<-read.csv("%s", sep="%s")''' % (os.path.join(gosetting.UPLOAD_FOLDER, parameters['universe']), parameters['duniverse'])
    robjects.r(universe_list)
    gene_list = '''genes<-read.csv("%s", sep="%s")''' % (os.path.join(gosetting.UPLOAD_FOLDER, parameters['study']), parameters['dstudy'])
    robjects.r(gene_list)
    gf = '''goFrame <- GOFrame(BostGO, organism = "%s")''' % parameters['organism']
    robjects.r(gf)
    gaf = '''goAllFrame <- GOAllFrame(goFrame)'''
    lg = robjects.r(gaf)
    #print(lg)
    gscf = '''gsc<-GeneSetCollection(goAllFrame,setType=GOCollection())'''
    lgc = robjects.r(gscf)
    #print(lgc)
    filename = parameters['association'].replace(" ", "_").replace(".csv",".txt")
    with open(os.path.join(gosetting.RESULTS_FOLDER, 'GO_'+filename), 'wt', newline='', encoding='utf-8') as output:
        fieldnames = ['GOID','Ontology','Pvalue','OddsRatio','ExpCount','Count','Size','Term','GeneIDs']
        output_write = csv.DictWriter(output, fieldnames=fieldnames, dialect='excel', delimiter='\t')
        output_write.writeheader()
        for i in parameters['ontology']:
            #print(i)
            para = '''params<-GSEAGOHyperGParams(name="My params",geneSetCollection=gsc,geneIds=genes,universeGeneIds=universe,pvalueCutoff=%f,testDirection="over",ontology="%s",conditional=FALSE)''' % (parameters['pvalue'], i)
            #print(para)
            robjects.r(para)
            proc = '''Over<-hyperGTest(params)'''
            robjects.r(proc)
            go_genes = '''geneIdsByCategory(Over)'''
            g = robjects.r(go_genes)
            g_dict = recurList(g)
            summary_res = '''head(summary(Over))'''
            s = robjects.r(summary_res)
            
            s_dict = recurList(s)
            result_dict = dict()
            for i0, i1, i2, i3, i4, i5, i6 in zip(s_dict['GO'+i+'ID'], s_dict['Pvalue'], s_dict['OddsRatio'], s_dict['ExpCount'], s_dict['Count'], s_dict['Size'], s_dict['Term']):
                geneids = ''
                if i0 in g_dict:
                    geneids_list = g_dict[i0]
                    geneids = ';'.join(geneids_list)
                result_dict = dict({'GOID': i0, 'Ontology': i, 'Pvalue': i1, 'OddsRatio':i2, 'ExpCount':i3, 'Count': i4, 'Size': i5, 'Term': i6, 'GeneIDs': geneids})
                output_write.writerow(result_dict)
        return ['GO_'+filename]
