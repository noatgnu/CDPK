from lxml import etree
#from xml.etree import ElementTree as etree
import os
import gzip
import csv
import re
from . import settingup
import MySQLdb
def allowed_file(filename):
    split_fn = filename.split('.')[-1]
    if split_fn:
       if split_fn in settingup.ALLOWED_EXTENSIONS:
           return 'True'
       else:
           return 'False'
    return 'False'

def check_column(input_path):
    condition = 0
    accession_col_name = ''
    with open(input_path, 'rt') as csvfile:
        reader = csv.DictReader(csvfile, dialect='excel', delimiter='\t')
        for row in reader:
            for k in row:
                if re.search('\|[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}\|', str(row[k])):
                    accession_col_name = k
                    condition = 1
            if condition == 0:
                if re.search('^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', str(row[k])):
                    accession_col_name = k
                    condition = 1
            break
    return accession_col_name
    
def uniprot_parser(parameters):
    list_entry = set()
    
    colm = check_column(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']))
    with open(os.path.join(settingup.UPLOAD_FOLDER, parameters['input']), 'rt') as infile:
        reader = csv.DictReader(infile, dialect='excel', delimiter='\t')
        for row in reader:
            result = re.search('([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', str(row[colm]))
            list_entry.add(result.group(0))
    
    ns = '{http://uniprot.org/uniprot}'
    sub_list = ''
    database_path = ''
    if 'database' in parameters:
        db = settingup.DB[parameters['database'][0].decode('utf-8')]
        database_path = os.path.join(settingup.APP_ROOT, db)
    go_found = False
    found_entry = False
    match_found = False
    relative_pos_acc = False
    relative_pos_name = False
    sub_loc = False
    recomname_loc = False
    goid_loc = False
    feature_loc = False
    feature_loc_seq = False
    dbr_id = False
    sub_listd = False
    go_idd = False
    go_termd = False
    seq_fd = False
    ns = '{http://uniprot.org/uniprot}'
    output = dict()
    protein = ''
    recommendedn = ''
    subcellular_loc = ''
    dbreference = ''
    feature_inf = ''
    position_seq = ''
    fieldname = ['accession']
    if 'entryname' in parameters:
        fieldname.append('entry name')
    if 'proteinname' in parameters:
        fieldname.append('protein name')
    if 'subcellularlocation' in parameters:
        fieldname.append('subcellular location')
    if 'goid' in parameters:
        fieldname.append('geneontology IDs')
    if 'goterm' in parameters:
        fieldname.append('geneontology terms')
    if 'sequence' in parameters:
        fieldname.append('sequence')
    if 'sequencefeatures' in parameters:
        fieldname.append('sequence features')
    if 'glycosylation' in parameters:
        fieldname.append('glycosylation sites')
    for k in fieldname:
        row[k] = ''
    with open(os.path.join(settingup.RESULTS_FOLDER, 'matched_'+parameters['input']), 'wt') as outfile:
        #fieldname = ['accession','entry name','protein name','subcellular location']
        output['Main Output'] = parameters['input']
        owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        owriter.writeheader()
        row = dict()   
        
        goid_list = ''
        goterm_list = ''
        seq_features = ''
        seq_position = ''
        ind_feature = ''
        gs_file = ''
        gswriter = ''
        gsrow = dict()
        goframe = ''
        gs_e = dict()
        
        if 'gsstats' in parameters:
            gs_file = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_'+parameters['input']), 'wt')
            gs_empty = open(os.path.join(settingup.RESULTS_FOLDER, 'gs_empty_'+parameters['input']), 'wt')
            
            
            gswriter = csv.DictWriter(gs_file, fieldnames=['Entry','Gocode','Geneontology IDs'], dialect='excel', delimiter='\t')
            gswriter.writeheader()
            gs_e_writer = csv.DictWriter(gs_empty, fieldnames=['Entry','Gocode','Geneontology IDs'], dialect='excel', delimiter='\t')
            gs_e_writer.writeheader()
                    
        for event, elem in etree.iterparse(gzip.open(database_path, 'rb'), events = ("start","end")):
            #print(event, elem.tag)
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and elem.text in list_entry and found_entry == True and event == "end":
                match_found = True
                relative_pos_name = True
                row['accession'] = elem.text
                print(elem.text)
                
                
            if elem.tag == ns+'name' and found_entry == True and match_found == True and relative_pos_name == True and event == "end" and 'entryname' in parameters:
                row['entry name'] = elem.text
                relative_pos_name = False
                
            if elem.tag == ns+'recommendedName' and match_found == True and found_entry == True and event == "start" and 'proteinname' in parameters:
                recomname_loc = True
                
            if elem.tag == ns+'fullName' and match_found == True and recomname_loc == True and event == "end" and 'proteinname' in parameters:
                row['protein name'] = elem.text
                
            if elem.tag == ns+'ecNumber' and match_found == True and recomname_loc == True and event == "end" and 'proteinname' in parameters:
                row['protein name'] += ' (EC '+ elem.text + ')'
                
            if elem.tag == ns+'recommendedName' and match_found == True and found_entry == True and event == "end" and 'proteinname' in parameters:
                recomname_loc = False
                
            if elem.tag == ns+'comment' and match_found == True and elem.get('type') == 'subcellular location' and found_entry == True and event == "start" and 'subcellularlocation' in parameters:
                sub_loc = True
                #print(etree.tostring(elem,pretty_print=True))
                
            if elem.tag == ns+'location' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            if elem.tag == ns+'topology' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            if elem.tag == ns+'text' and match_found == True and found_entry == True and sub_loc == True and event == "end" and 'subcellularlocation' in parameters:
                if sub_listd == True:
                    sub_list = sub_list+';'
                sub_list = sub_list+elem.text
                sub_listd = True
            
            if elem.tag == ns+'comment' and match_found == True and elem.get('type') == 'subcellular location' and found_entry == True and event == "end" and 'subcellularlocation' in parameters:
                sub_loc = False
                row['subcellular location'] = sub_list
                sub_list = ''
                sub_listd = False
                #print(etree.tostring(elem, pretty_print=True))
            if elem.tag == ns+'dbReference' and match_found == True and found_entry == True and elem.get('type') == 'GO' and event == "start":
                go_found = True
                #print(go_idd)
                if 'goid' in parameters or 'gsstats' in parameters:
                    goid_loc = True
                    if go_idd == True:
                        goid_list = goid_list + ';' + elem.get('id')
                    else:
                        goid_list = elem.get('id')
                        go_idd = True
                     
                if 'goterm' in parameters:
                    goframe = elem.get('id')
                
            if elem.tag == ns+'property' and match_found == True and found_entry == True and goid_loc == True and 'goterm' in parameters:
                if elem.get('type') == 'term':
                    if go_termd == True:
                        goterm_list = goterm_list+';'
                        
                    goterm_list = goterm_list+goframe+'('+elem.get('value')+')'
                    go_termd = True
            if elem.tag == ns+'dbReference' and match_found == True and found_entry == True and elem.get('type') == 'GO' and event == "end":
                if 'goid' in parameters or 'gsstats' in parameters:
                    goid_loc = False
                
            if elem.tag == ns+'feature' and match_found == True and found_entry == True and event == "start":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    feature_loc = True
                    if seq_fd == True:
                        seq_features = seq_features + ';'
                    seq_fd = True
                    if 'description' in elem.attrib:
                        ind_feature = elem.get('type')+ ': ' + elem.get('description')
                    else:
                        ind_feature = elem.get('type')
                    
                    
            if elem.tag == ns+'location' and match_found == True and found_entry == True and event == "start":
                if 'goid' in parameters or 'gsstats' in parameters:
                    feature_loc_seq = True
                
            if elem.tag == ns+'begin' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end" and 'sequencefeatures' in parameters:
                
                if not 'position' in elem.attrib:
                    seq_position = 'unknown'
                        
                else:
                    seq_position = elem.get('position')
                                
            if elem.tag == ns+'end' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end" and 'sequencefeatures' in parameters:
                
                if not 'position' in elem.attrib:
                    seq_position = seq_position +'-'+ 'unknown'
                        
                else:
                    seq_position = seq_position +'-'+elem.get('position')
                
            if elem.tag == ns+'position' and match_found == True and found_entry == True and feature_loc_seq == True and event == "end":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    if not 'position' in elem.attrib:
                        seq_position = 'unknown'
                            
                    else:
                        seq_position = elem.get('position')
                    
            if elem.tag == ns+'location' and match_found == True and found_entry == True and event == "end":
                if 'sequencefeatures' in parameters or 'glycosylation' in parameters:
                    feature_loc_seq = False
                
            if elem.tag == ns+'feature' and match_found == True and found_entry == True and event == "end":
                if 'goid' in parameters or 'gsstats' in parameters:
                    feature_loc = False
                    seq_features = seq_features + ind_feature+'('+seq_position+')'
                
            if elem.tag == ns+'sequence' and match_found == True and found_entry == True and event == "end" and 'sequence' in parameters:
                row['sequence'] = elem.text
                
            if elem.tag == ns+'entry' and match_found == True and found_entry == True and event == "end":
                
                #print(goid_list)
                if go_found == True:
                    output['GSSTATS'] = 'gs_'+parameters['input']
                    if 'goid' in parameters:
                        row['geneontology IDs'] = goid_list
                    if 'goterm' in parameters:
                        row['geneontology terms'] = goterm_list
                else:
                    output['GSSTATS Entry without any GOIDs'] = 'gs_empty_'+parameters['input']
                    if 'goid' in parameters:
                        row['geneontology IDs'] = ''
                    if 'goterm' in parameters:
                        row['geneontology terms'] = ''
                    gs_e['Entry'] = row['accession']
                    gs_e['Gocode'] = 'IEA'
                    gs_e['Geneontology IDs'] = ''
                    gs_e_writer.writerow(gs_e)
                    
                    
                if 'sequencefeatures' in parameters:
                    row['sequence features'] = seq_features
                if 'glycosylation' in parameters:
                    row['glycosylation sites'] = seq_features
                if 'gsstats' in parameters:
                    if  'entry name' in row:
                        if not row['entry name'] == '':
                            gsrow['Entry'] = row['entry name']
                        else:
                            gsrow['Entry'] = row['accession']
                    else:
                        gsrow['Entry'] = row['accession']
                    goid_ar = goid_list.split(';')
                    #print(goid_ar)
                    for goid in goid_ar:
                        gsrow['Geneontology IDs'] = goid
                        gsrow['Gocode'] = 'IEA'
                        #print(gsrow)
                        gswriter.writerow(gsrow)
                owriter.writerow(row)
                list_entry.remove(row['accession'])    
                print(len(list_entry))
                goid_list = ''
                goterm_list = ''
                seq_features = ''
                go_idd = False
                go_termd = False
                seq_fd = False
                gs_e = dict()
                gsrow = dict()
                elem.clear()
                go_found = False
                found_entry = False
                match_found = False
                for k in fieldname:
                    row[k] = ''
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
                if len(list_entry) == 0:
                    if 'gsstats' in parameters:
                        gs_file.close()
                        gs_empty.close()
                    return output                    
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        if 'gsstats' in parameters:
            gs_file.close()
            gs_empty.close()
        if not len(list_entry) == 0:
            output['No Match'] = 'no_match_'+parameters['input']
            field_mini = ['accession']
            with open(os.path.join(settingup.RESULTS_FOLDER, 'no_match_'+parameters['input']), 'wt') as gs_file_no_hit:
                no_match = csv.DictWriter(gs_file_no_hit, fieldnames=field_mini, dialect='excel', delimiter='\t')
                no_match.writeheader()
                no_match_row = dict()
                for i in list_entry:
                    no_match_row['accession'] = i
                    no_match.writerow(no_match_row)
    return output
    
def compdb(filepath='uniprot_sprot.xml.gz', output='compact_uniprot_sprot.xml.gz'):
    
    
    with gzip.open(output, 'wb', 9) as outfile:
        found_entry = False
        relative_pos_acc = False
        relative_pos_name = False
        sub_loc = False
        recomname_loc = False
        goid_loc = False
        feature_loc = False
        feature_loc_seq = False
        dbr_id = False
        ns = '{http://uniprot.org/uniprot}'
        entry = etree.Element('entry')
        protein = ''
        recommendedn = ''
        subcellular_loc = ''
        #subr_loc = etree.Element('subcellularlocation')
        dbreference = ''
        feature_inf = ''
        position_seq = ''
        found_organism = False
        #fieldname = ['accession','entry name','protein name','subcellular location']
        #owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        #owriter.writeheader()
        outfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        outfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        #row = dict()    
        for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
            #print(event, etree.tostring(elem))
            #if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and event == 'start':
                #print(etree.tostring(elem))
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and found_entry == True and event == "end":
                
                accession = etree.SubElement(entry, 'accession')
                accession.text = elem.text
                #print(elem.text)
                #print(etree.tostring(entry))
                relative_pos_name = True
            if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
                e_name = etree.SubElement(entry,'name')
                e_name.text = elem.text
                
                relative_pos_name = False
            if elem.tag == ns+'protein' and found_entry == True and event == "start":
                protein = etree.SubElement(entry,'protein')
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
                recomname_loc = True
                recommendedn = etree.SubElement(protein, 'recommendedName')
            if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
                fullname = etree.SubElement(recommendedn, 'fullName')
                fullname.text = elem.text
            if elem.tag == ns+'ecNumber' and recomname_loc == True and event == "end":
                ecnm = etree.SubElement(recommendedn, 'ecNumber')
                ecnm.text = elem.text
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
                recomname_loc = False
            if elem.tag == ns+'organism' and found_entry == True and event == "start":
                found_organism = True
                organism = etree.SubElement(entry, 'organism')
            if elem.tag == ns+'name' and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
                org_name = etree.SubElement(organism, 'name', {'type':'scientific'})
                org_name.text = elem.text
            if elem.tag == ns+'organism' and found_entry == True and event == "end":
                found_organism = False
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
                sub_loc = True
                subcellular_loc = etree.SubElement(entry, 'comment', {'type': 'subcellular location'})
            if elem.tag == ns+'subcellularLocation' and found_entry == True and event == "start":
                
                subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
                  
            if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
                location = etree.SubElement(subr_loc, 'location')
                location.text = elem.text
                #print(elem.text)
            if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
                topology = etree.SubElement(subr_loc, 'topology')
                topology.text = elem.text
            if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
                
                topology = etree.SubElement(subr_loc, 'orientation')
                topology.text = elem.text
                
            if elem.tag == ns+'text' and found_entry == True and sub_loc == True and event == "end":
                text = etree.SubElement(subcellular_loc, 'text')
                text.text = elem.text
                
                #print(elem.text)
            #if elem.tag == ns+'subcellular location' and found_entry == True and event == "end":
                #subcellular_loc.append(subr_loc)
                #subr_loc = etree.Element('subcellular location')
                
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
                sub_loc = False
                
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
                goid_loc = True
                dbre_at = {'type':elem.get('type'), 'id':elem.get('id')}
                dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
                
            if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
                dbre_at = {'type':elem.get('type'), 'value':elem.get('value')}
                dbrprop = etree.SubElement(dbReference, 'property', dbre_at)
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
                goid_loc = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "start":
                feature_loc = True
                fattribute = {'type':elem.get('type')}
                if elem.get('description'):
                    fattribute['description'] = elem.get('description')
                    
                feature_inf = etree.SubElement(entry, 'feature', fattribute)
                
                
            if elem.tag == ns+'location' and found_entry == True and event == "start":
                position_seq = etree.SubElement(feature_inf, 'location')
                feature_loc_seq = True
            if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
                b_start = dict()
                if not 'position' in elem.attrib:
                    b_start['status'] = 'unknown'
                        
                else:
                    b_start['position'] = elem.get('position')
                
                b_pos = etree.SubElement(position_seq, 'begin', b_start)
                
            if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
                b_end = dict()
                if not 'position' in elem.attrib:
                    b_end['status'] = 'unknown'
                        
                else:
                    b_end['position'] = elem.get('position')
                e_pos = etree.SubElement(position_seq, 'end', b_end)
                
            if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
                b_direct = dict()
                if not 'position' in elem.attrib:
                    b_direct['status'] = 'unknown'
                        
                else:
                    b_direct['position'] = elem.get('position')
                pos = etree.SubElement(position_seq, 'position', b_direct)
                
            if elem.tag == ns+'location' and found_entry == True and event == "end":
                feature_loc_seq = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "end":
                feature_loc = False
            if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
                sattribute = dict()
                if elem.get('length'):
                    sattribute['length'] = elem.get('length')
                if elem.get('mass'):
                    sattribute['mass'] = elem.get('mass')
                sequence = etree.SubElement(entry, 'sequence', sattribute)
                
                sequence.text = elem.text
            if elem.tag == ns+'entry' and found_entry == True and event == "end":
                found_entry = False
                
                #print(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                outfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                elem.clear()
                entry = etree.Element('entry')
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        outfile.write(b'</uniprot>')

def compdb_sql(filepath='uniprot_sprot.xml.gz', dtype='sp'):
    count = 0
    p1n = 'compactUniprot'
    ac1 = p1n+dtype+'ACC'
    fn1 = p1n+dtype+'FullName'
    ec1 = p1n+dtype+'EC'
    or1 = p1n+dtype+'Organism'
    sl1 = p1n+dtype+'SubLoc'
    go1 = p1n+dtype+'GO'
    ft1 = p1n+dtype+'FT'
    se1 = p1n+dtype+'Seq'
    dbuser = settingup.DATABASES['default']['USER']
    dbpass = str(settingup.DATABASES['default']['PASSWORD'])
    dbhost = settingup.DATABASES['default']['HOST']
    db = settingup.DATABASES['default']['DB']
    dbcon = MySQLdb.connect(user=dbuser, passwd=dbpass, host=dbhost, db=db)
    cursor = dbcon.cursor()
    dacc = """DROP TABLE IF EXISTS %s;""" % ac1
    cursor.execute(dacc)
    dfn = """DROP TABLE IF EXISTS %s;""" % fn1
    cursor.execute(dfn)
    dec = """DROP TABLE IF EXISTS %s;""" % ec1
    cursor.execute(dec)
    dorg = """DROP TABLE IF EXISTS %s;""" % or1
    cursor.execute(dorg)
    dsubl = """DROP TABLE IF EXISTS %s;""" % sl1
    cursor.execute(dsubl)
    dgo = """DROP TABLE IF EXISTS %s;""" % go1
    cursor.execute(dgo)
    dft = """DROP TABLE IF EXISTS %s;""" % ft1
    cursor.execute(dft)
    dseq = """DROP TABLE IF EXISTS %s;""" % se1
    cursor.execute(dseq)
    tbc_acc = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, accession char(20) NOT NULL);""" % ac1
    cursor.execute(tbc_acc)
    tbc_fn = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, fullName text NOT NULL);""" % fn1
    cursor.execute(tbc_fn)
    tbc_ec = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, EC char(20) NOT NULL);""" % ec1
    cursor.execute(tbc_ec)
    tbc_org = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, organismName text NOT NULL);""" % or1
    cursor.execute(tbc_org)
    tbc_sub_loc = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, subType char(20) NOT NULL, location text NOT NULL);""" % sl1
    cursor.execute(tbc_sub_loc)
    tbc_go = """CREATE TABLE IF NOT EXISTS %s (entryNumber int(20) NOT NULL, entryName char(20) NOT NULL, GOID char(20) NOT NULL, evidence char(255) NOT NULL, GOTerm text NOT NULL);""" % go1
    cursor.execute(tbc_go)
    tbc_ft = """CREATE TABLE IF NOT EXISTS %s (entryNumber int(20) NOT NULL, entryName char(20) NOT NULL, FTType text NOT NULL, position char(50) NOT NULL, description text NOT NULL);""" % ft1
    cursor.execute(tbc_ft)
    tbc_seq = """CREATE TABLE IF NOT EXISTS %s (entryName char(20) NOT NULL, sequence text NOT NULL);""" % se1
    cursor.execute(tbc_seq)
    #sln = 0
    gon = 0
    ftn = 0
    
    ins_p = """INSERT INTO """
    up_p = """UPDATE """
    
    up_go = """ SET evidence = %s WHERE entryNumber = """
    up_go2 = """ SET GOTerm = %s WHERE entryNumber = """
    up_ft = """ SET position = %s WHERE entryNumber = """
    up_ft2 = """ SET description = %s WHERE entryNumber = """
    ins_acc = """ (entryName, accession) VALUES (%s, %s);"""
    ins_fn = """ (entryName, fullName) VALUES (%s, %s);"""
    ins_ec = """ (entryName, EC) VALUES (%s, %s);"""
    ins_org = """ (entryName, organismName) VALUES (%s, %s);"""
    ins_sub_loc = """ (entryName, subType, location) VALUES (%s, %s, %s);"""
    ins_go = """ (entryNumber, entryName, GOID, evidence, GOTerm) VALUES (%s, %s, %s, %s, %s);"""
    ins_ft = """ (entryNumber, entryName, FTType, position, description) VALUES (%s, %s, %s, %s, %s);"""
    ins_seq = """ (entryName, sequence) VALUES (%s, %s);"""
    #with gzip.open(output, 'wb', 9) as outfile:
    found_entry = False
    relative_pos_acc = False
    relative_pos_name = False
    sub_loc = False
    recomname_loc = False
    goid_loc = False
    feature_loc = False
    feature_loc_seq = False
    dbr_id = False
    found_organism = False
    ns = '{http://uniprot.org/uniprot}'
    #entry = etree.Element('entry')
    protein = ''
    recommendedn = ''
    subcellular_loc = ''
    #subr_loc = etree.Element('subcellularlocation')
    dbreference = ''
    feature_inf = ''
    position_seq = ''
    access_entry = ''
    accession_list = []
    #sl_stuff = ''
    go_stuff = ''
    
    posf = ''
    #fieldname = ['accession','entry name','protein name','subcellular location']
    #owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
    #owriter.writeheader()
    #outfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
    #outfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
    #row = dict()    
    for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
        #print(event, etree.tostring(elem))
        #if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and event == 'start':
            #print(etree.tostring(elem))
        
        if elem.tag == ns+'entry' and event == "start":
            found_entry = True
            accession_list = []
            
        if elem.tag == ns+'accession' and found_entry == True and event == "end":
            
            #accession = etree.SubElement(entry, 'accession')
            #accession.text = elem.text
            accession_list.append(elem.text)
            
            #print(elem.text)
            #print(etree.tostring(entry))
            relative_pos_name = True
        if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
            #e_name = etree.SubElement(entry,'name')
            #e_name.text = elem.text
            access_entry = elem.text
            for acc in accession_list:
                cursor.execute(ins_p+ac1+ins_acc, (access_entry, acc,))
                dbcon.commit()
            relative_pos_name = False
        #if elem.tag == ns+'protein' and found_entry == True and event == "start":
            #protein = etree.SubElement(entry,'protein')
        if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
            recomname_loc = True
            #recommendedn = etree.SubElement(protein, 'recommendedName')
        if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
            #fullname = etree.SubElement(recommendedn, 'fullName')
            #fullname.text = elem.text
            cursor.execute(ins_p+fn1+ins_fn, (access_entry, elem.text))
            dbcon.commit()
        if elem.tag == ns+'ecNumber' and recomname_loc == True and event == "end":
            #ecnm = etree.SubElement(recommendedn, 'ecNumber')
            #ecnm.text = elem.text
            cursor.execute(ins_p+ec1+ins_ec, (access_entry, elem.text))
            dbcon.commit()
        if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
            recomname_loc = False
        if elem.tag == ns+'organism' and found_entry == True and event == "start":
            found_organism = True
            #organism = etree.SubElement(entry, 'organism')
        if elem.tag == ns+'name' and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
            #org_name = etree.SubElement(organism, 'name', {'type':'scientific'})
            #org_name.text = elem.text
            cursor.execute(ins_p+or1+ins_org, (access_entry, elem.text))
            dbcon.commit()
            found_organism = False
        if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
            sub_loc = True
            #subcellular_loc = etree.SubElement(entry, 'comment', {'type': 'subcellular location'})
        #if elem.tag == ns+'subcellularLocation' and found_entry == True and event == "start":
            
            #subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
              
        if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
            #print(elem.text)
            #location = etree.SubElement(subr_loc, 'location')
            #location.text = elem.text
            
            cursor.execute(ins_p+sl1+ins_sub_loc, (access_entry, 'location', elem.text))
            dbcon.commit()
            
        if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
            #print(elem.text)
            #topology = etree.SubElement(subr_loc, 'topology')
            #topology.text = elem.text
            
            cursor.execute(ins_p+sl1+ins_sub_loc, (access_entry, 'topology', elem.text))
            dbcon.commit()
            
        if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
            #print(elem.text)
            #topology = etree.SubElement(subr_loc, 'orientation')
            #topology.text = elem.text
            
            cursor.execute(ins_p+sl1+ins_sub_loc, (access_entry, 'location', elem.text))
            dbcon.commit()
        #if elem.tag == ns+'text' and found_entry == True and sub_loc == True and event == "end":
            #text = etree.SubElement(subcellular_loc, 'text')
            #text.text = elem.text
            
            #print(elem.text)
        #if elem.tag == ns+'subcellular location' and found_entry == True and event == "end":
            #subcellular_loc.append(subr_loc)
            #subr_loc = etree.Element('subcellular location')
            
        if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
            sub_loc = False
            
        if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
            gon += 1
            goid_loc = True
            #dbre_at = {'type':elem.get('type'), 'id':elem.get('id')}
            #dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
            
            go_stuff = gon
            cursor.execute(ins_p+go1+ins_go, (str(go_stuff), access_entry, elem.get('id'), '', ''))
            dbcon.commit()
        if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
            #dbre_at = {'type':elem.get('type'), 'value':elem.get('value')}
            #dbrprop = etree.SubElement(dbReference, 'property', dbre_at)
            if elem.get('type') == 'term':
                cursor.execute(up_p+go1+up_go+str(go_stuff)+';', (elem.get('value'),))
                dbcon.commit()
            if elem.get('type') == 'evidence':
                #print(elem.get('value'))
                cursor.execute(up_p+go1+up_go2+str(go_stuff)+';', (elem.get('value'),))
                dbcon.commit()
        if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
            goid_loc = False
            
        if elem.tag == ns+'feature' and found_entry == True and event == "start":
            ftn += 1
            psof = ''
            feature_loc = True
            #fattribute = dict()
            #fattribute['type'] = elem.get('type')
            ft_stuff = ftn
            cursor.execute(ins_p+ft1+ins_ft, (str(ft_stuff), access_entry, elem.get('type'), '', ''))
            dbcon.commit()
            if 'description' in elem.attrib:
                #fattribute['description'] = elem.get('description')
                
                cursor.execute(up_p+ft1+up_ft2+str(ft_stuff)+';', (elem.get('description'),))
                dbcon.commit()
            #feature_inf = etree.SubElement(entry, 'feature', fattribute)
            
            
        if elem.tag == ns+'location' and found_entry == True and event == "start":
            #position_seq = etree.SubElement(feature_inf, 'location')
            feature_loc_seq = True
            posf = ''
        if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
            #b_start = dict()
            if not 'position' in elem.attrib:
                #b_start['status'] = 'unknown'
                posf = posf + 'unknown'
                    
            else:
                #b_start['position'] = elem.get('position')
                posf = posf + elem.get('position')
            #b_pos = etree.SubElement(position_seq, 'begin', b_start)
            
        if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
            #b_end = dict()
            if not 'position' in elem.attrib:
                #b_end['status'] = 'unknown'
                posf = posf + '-unknown'
            else:
                #b_end['position'] = elem.get('position')
                posf = posf +'-'+ elem.get('position')
            cursor.execute(up_p+ft1+up_ft+str(ft_stuff)+';', (posf,))
            dbcon.commit()
            #e_pos = etree.SubElement(position_seq, 'end', b_end)
            
        if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
            #b_direct = dict()
            if not 'position' in elem.attrib:
                #b_direct['status'] = 'unknown'
                cursor.execute(up_p+ft1+up_ft+str(ft_stuff)+';', ('unknown',))
                dbcon.commit()
            else:
                #b_direct['position'] = elem.get('position')
                cursor.execute(up_p+ft1+up_ft+str(ft_stuff)+';', (elem.get('position'),))
                dbcon.commit()
            #pos = etree.SubElement(position_seq, 'position', b_direct)
            
        if elem.tag == ns+'location' and found_entry == True and event == "end":
            feature_loc_seq = False
            
        if elem.tag == ns+'feature' and found_entry == True and event == "end":
            feature_loc = False
            
        if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
            #sattribute = dict()
            #if elem.get('length'):
                #sattribute['length'] = elem.get('length')
            #if elem.get('mass'):
                #sattribute['mass'] = elem.get('mass')
            #sequence = etree.SubElement(entry, 'sequence', sattribute)
            #sequence.text = elem.text
            cursor.execute(ins_p+se1+ins_seq, (access_entry, elem.text))
            dbcon.commit()
        if elem.tag == ns+'entry' and found_entry == True and event == "end":
            count+=1
            
            print(count)
            found_entry = False
            access_entry = ''
            #print(etree.tostring(entry, encoding='utf-8', pretty_print=True))
            #outfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
            elem.clear()
            #entry = etree.Element('entry')
            while elem.getprevious() is not None:
                del elem.getparent()[0]
        if elem.tag == ns+'entry' and found_entry == False and event == "end":
            found_entry = False
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
                
    #outfile.write(b'</uniprot>')
    cursor.close()
    dbcon.close()
 
def compdb_taxa(filepath='uniprot_sprot.xml.gz', output='compact_uniprot_sprot.xml.gz', taxa_list=['Vertebrata', 'Mammalia', 'Fungi', 'Eukaryota']):
    taxa_dict = dict()
    for i in taxa_list:
        taxa_dict[i] = False
    count = 0
    with gzip.open(output, 'wb', 9) as outfile, gzip.open('vertebrata_'+output, 'wb', 9) as vboutfile, gzip.open('mammalia_'+output, 'wb', 9) as mmoutfile, gzip.open('fungi_'+output, 'wb', 9) as fgoutfile, gzip.open('eukaryota_'+output, 'wb', 9) as ekoutfile:
        found_entry = False
        relative_pos_acc = False
        relative_pos_name = False
        sub_loc = False
        recomname_loc = False
        goid_loc = False
        feature_loc = False
        feature_loc_seq = False
        dbr_id = False
        ns = '{http://uniprot.org/uniprot}'
        entry = etree.Element('entry')
        protein = ''
        recommendedn = ''
        subcellular_loc = ''
        #subr_loc = etree.Element('subcellularlocation')
        dbreference = ''
        feature_inf = ''
        position_seq = ''
        found_organism = False
        #fieldname = ['accession','entry name','protein name','subcellular location']
        #owriter = csv.DictWriter(outfile, fieldnames=fieldname, dialect='excel', delimiter='\t')
        #owriter.writeheader()
        outfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        outfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        #row = dict()
        
        vboutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        vboutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
            
        mmoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        mmoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        
        fgoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        fgoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        
        ekoutfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        ekoutfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
            
        for event, elem in etree.iterparse(gzip.open(filepath, 'rb'), events = ("start","end")):
            #print(event, etree.tostring(elem))
            #if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and event == 'start':
                #print(etree.tostring(elem))
            
            if elem.tag == ns+'entry' and event == "start":
                found_entry = True 
                
            if elem.tag == ns+'accession' and found_entry == True and event == "end":
                
                accession = etree.SubElement(entry, 'accession')
                accession.text = elem.text
                #print(elem.text)
                #print(etree.tostring(entry))
                relative_pos_name = True
            if elem.tag == ns+'name' and found_entry == True and relative_pos_name == True and event == "end":
                e_name = etree.SubElement(entry,'name')
                e_name.text = elem.text
                
                relative_pos_name = False
            if elem.tag == ns+'protein' and found_entry == True and event == "start":
                protein = etree.SubElement(entry,'protein')
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "start":
                recomname_loc = True
                recommendedn = etree.SubElement(protein, 'recommendedName')
            if elem.tag == ns+'fullName' and recomname_loc == True and event == "end":
                fullname = etree.SubElement(recommendedn, 'fullName')
                fullname.text = elem.text
            if elem.tag == ns+'ecNumber' and recomname_loc == True and event == "end":
                ecnm = etree.SubElement(recommendedn, 'ecNumber')
                ecnm.text = elem.text
            if elem.tag == ns+'recommendedName' and found_entry == True and event == "end":
                recomname_loc = False
            if elem.tag == ns+'organism' and found_entry == True and event == "start":
                found_organism = True
                organism = etree.SubElement(entry, 'organism')
            if elem.tag == ns+'name' and found_entry == True and elem.get('type') == 'scientific' and found_organism == True and event == "end":
                org_name = etree.SubElement(organism, 'name', {'type':'scientific'})
                org_name.text = elem.text
            if elem.tag == ns+'taxon' and found_organism == True and found_entry == True and event == "start" and elem.text in taxa_dict:
                taxa_dict[elem.text] = True
            if elem.tag == ns+'organism' and found_entry == True and event == "end":
                found_organism = False
            
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "start":
                sub_loc = True
                subcellular_loc = etree.SubElement(entry, 'comment', {'type': 'subcellular location'})
            
            if elem.tag == ns+'subcellularLocation' and found_entry == True and event == "start":
                
                subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
                  
            if elem.tag == ns+'location' and found_entry == True and sub_loc == True and event == "end":
                location = etree.SubElement(subr_loc, 'location')
                location.text = elem.text
                #print(elem.text)
            if elem.tag == ns+'topology' and found_entry == True and sub_loc == True and event == "end":
                topology = etree.SubElement(subr_loc, 'topology')
                topology.text = elem.text
            if elem.tag == ns+'orientation' and found_entry == True and sub_loc == True and event == "end":
                
                topology = etree.SubElement(subr_loc, 'orientation')
                topology.text = elem.text
                
            if elem.tag == ns+'text' and found_entry == True and sub_loc == True and event == "end":
                text = etree.SubElement(subcellular_loc, 'text')
                text.text = elem.text
                
                #print(elem.text)
            #if elem.tag == ns+'subcellular location' and found_entry == True and event == "end":
                #subcellular_loc.append(subr_loc)
                #subr_loc = etree.Element('subcellular location')
                
            if elem.tag == ns+'comment' and elem.get('type') == 'subcellular location' and found_entry == True and event == "end":
                sub_loc = False
                
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "start":
                goid_loc = True
                dbre_at = {'type':elem.get('type'), 'id':elem.get('id')}
                dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
                
            if elem.tag == ns+'property' and found_entry == True and goid_loc == True and event == "end":
                dbre_at = {'type':elem.get('type'), 'value':elem.get('value')}
                dbrprop = etree.SubElement(dbReference, 'property', dbre_at)
            if elem.tag == ns+'dbReference' and found_entry == True and elem.get('type') == 'GO' and event == "end":
                goid_loc = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "start":
                feature_loc = True
                fattribute = {'type':elem.get('type')}
                if elem.get('description'):
                    fattribute['description'] = elem.get('description')
                    
                feature_inf = etree.SubElement(entry, 'feature', fattribute)
                
                
            if elem.tag == ns+'location' and found_entry == True and event == "start":
                position_seq = etree.SubElement(feature_inf, 'location')
                feature_loc_seq = True
            if elem.tag == ns+'begin' and found_entry == True and feature_loc_seq == True and event == "end":
                b_start = dict()
                if not 'position' in elem.attrib:
                    b_start['status'] = 'unknown'
                        
                else:
                    b_start['position'] = elem.get('position')
                
                b_pos = etree.SubElement(position_seq, 'begin', b_start)
                
            if elem.tag == ns+'end' and found_entry == True and feature_loc_seq == True and event == "end":
                b_end = dict()
                if not 'position' in elem.attrib:
                    b_end['status'] = 'unknown'
                        
                else:
                    b_end['position'] = elem.get('position')
                e_pos = etree.SubElement(position_seq, 'end', b_end)
                
            if elem.tag == ns+'position' and found_entry == True and feature_loc_seq == True and event == "end":
                b_direct = dict()
                if not 'position' in elem.attrib:
                    b_direct['status'] = 'unknown'
                        
                else:
                    b_direct['position'] = elem.get('position')
                pos = etree.SubElement(position_seq, 'position', b_direct)
                
            if elem.tag == ns+'location' and found_entry == True and event == "end":
                feature_loc_seq = False
                
            if elem.tag == ns+'feature' and found_entry == True and event == "end":
                feature_loc = False
            if elem.tag == ns+'sequence' and found_entry == True and event == "end" and 'checksum' in elem.attrib:
                sattribute = dict()
                if elem.get('length'):
                    sattribute['length'] = elem.get('length')
                if elem.get('mass'):
                    sattribute['mass'] = elem.get('mass')
                sequence = etree.SubElement(entry, 'sequence', sattribute)
                
                sequence.text = elem.text
            if elem.tag == ns+'entry' and found_entry == True and event == "end":
                found_entry = False
                
                #print(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                outfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                if taxa_dict['Vertebrata'] == True:
                    vboutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Vertebrata'] = False
                if taxa_dict['Mammalia'] == True:
                    mmoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Mammalia'] = False
                if taxa_dict['Fungi'] == True:
                    fgoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Fungi'] = False
                if taxa_dict['Eukaryota'] == True:
                    ekoutfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
                    taxa_dict['Eukaryota'] = False
                count += 1
                print(count, end='\r')
                elem.clear()
                entry = etree.Element('entry')
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
            if elem.tag == ns+'entry' and found_entry == False and event == "end":
                found_entry = False
                elem.clear()
                while elem.getprevious() is not None:
                    del elem.getparent()[0]
        outfile.write(b'</uniprot>')
        vboutfile.write(b'</uniprot>')
        mmoutfile.write(b'</uniprot>')
        fgoutfile.write(b'</uniprot>')
        ekoutfile.write(b'</uniprot>')
