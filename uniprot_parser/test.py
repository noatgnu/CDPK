import gzip
from lxml import etree
import re
def uniprot_dat_compacter(fileinput='uniprot_sprot_fungi.dat.gz',fileoutput='uniprot_sprot_fungi.xml.gz'):
    with gzip.open(fileinput, 'rt') as inputfile, gzip.open(fileoutput, 'wb', 9) as outputfile:
        
        outputfile.write(b'<?xml version="1.0" encoding="UTF-8"?>')
        outputfile.write(b'<uniprot xmlns="http://uniprot.org/uniprot" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd">')
        count = 0
        entry_name = ''
        acc_list = []
        subcell_loc = []
        goid_list = []
        goterm_list = []
        fn = ''
        ecn = ''
        FT_array = []
        CC_string = ''
        seq_aa = ''
        seq_mw = ''
        start_entry = False
        seq_loc = False
        seq = ''
        subcc = False
        go_found = False
        ft_found = False
        
        for line in inputfile:
            
            if line.startswith('ID   '):
                entry_name = ''
                start_entry = True
                res = re.search('[^\s]*', line[5:])
                entry_name = res.group(0)
                #print(entry_name)
            if line.startswith('AC   ') and start_entry == True:
                for i in line[5:].split(';'):
                    if not i == '':
                        acc_list.append(i.lstrip().rstrip())
                        #print(acc_list)
            if line.startswith('DE   ') and start_entry == True:
                if line[5:].startswith('RecName: '):
                    fn_re = re.search('Full=([^;|\{]*)', line)
                    if res is not None:
                        fn = fn_re.group(1)
                ec = re.search('EC=([^;|\{]*)', line)
                #print(fn)
                if ec is not None:
                    ecn = ec.group(1)
                #print(ecn)
            if line.startswith('CC   ') and start_entry == True:
                
                CC_string = CC_string + line[5:].lstrip()
                CC_list = CC_string.split('-!-')
                for cc in CC_list:
                    c = cc.lstrip().rstrip()
                    if c.startswith('SUBCELLULAR LOCATION: '):
                        subcc = True
                        sub_loc = c.replace('SUBCELLULAR LOCATION: ', '').replace(';', ',').strip().replace('\n', '').replace('\r', '')
                        #print(sub_loc)
                        sub_1 = re.sub('-----.*-----', '', sub_loc)
                        #print(sub_1)
                        sub_list = re.sub('\{.*?\}', '', sub_1)
                        #print(sub_list)
                        subcell_loc = sub_list.split('.')
                        
                        #print(subcell_loc)
            if line.startswith('DR   ') and start_entry == True:
                if line[5:].startswith('GO;'):
                    go_found = True
                    goid_re = re.search('(GO:\w*);', line)
                    goid_list.append(goid_re.group(1))
                    #print(goid_list)
                    goterm_re = re.search('([C|F|P]:.*);', line)
                    goterm_list.append(goterm_re.group(1))
            if line.startswith('FT   ') and start_entry == True:
                ft_found = True
                
                if not line[5:].startswith(' '):
                    ft_frame = ''
                    ft_frame = line[5:].strip().replace('\n', '').replace('\r', '')
                if line[5:].startswith(' '):
                    ft_frame = ft_frame + ' ' + line[5:].strip().replace('\n', '').replace('\r', '')
                #print(ft_frame)
                FT_re = re.search('(^\w*)\s*(\d*)\s*(\d*)\s*(.+)', ft_frame)
                #print(FT_re.group(4))
                FTT_re = re.sub('\{.*?\}\.', '', FT_re.group(4))
                #print(FTT_re)
                FT4_re = re.sub('/FTId=.*?\.', '', FTT_re)
                FT_array.append([FT_re.group(1),FT_re.group(2),FT_re.group(3),FT4_re])
                #print(FT_array)
            if line.startswith('SQ   ') and start_entry == True:
                seq = ''
                seq_aa = ''
                seq_mw = ''
                if line[5:].startswith('SEQUENCE'):
                    seq_loc = True
                    seqinf_str = line[5:].replace('SEQUENCE', '')
                    aa = re.search('(\d*)\sAA', seqinf_str)
                    if aa is not None:
                        seq_aa = str(aa.group(1))
                    mw = re.search('(\d*)\sMW', seqinf_str)
                    if mw is not None:
                        seq_mw = str(mw.group(1))
            if line.startswith('     ') and seq_loc == True and start_entry == True:
                seq = seq+line[5:].strip().replace('\r', '').replace('\n', '').replace(' ','')
                #print(seq)
            if line.startswith('//') and start_entry == True:
                #print(line)
                
                entry = etree.Element('entry')
                for acc in acc_list:
                    if not acc == '':
                        accession = etree.SubElement(entry, 'accession')
                        accession.text = acc
                e_name = etree.SubElement(entry,'name')
                e_name.text = entry_name
                
                protein = etree.SubElement(entry,'protein')
                recommendedn = etree.SubElement(protein, 'recommendedName')
                fullname = etree.SubElement(recommendedn, 'fullName')
                fullname.text = fn
                if not ecn == '':
                    ecnm = etree.SubElement(recommendedn, 'ecNumber')
                    ecnm.text = ecn
                #print(etree.tostring(entry))
                if subcc == True:
                    subcellular_loc = etree.SubElement(entry, 'comment', {'type':'subcellular location'})
                    
                    for sc in subcell_loc:
                        if not sc == '':
                            if not sc.strip().startswith('Note='):
                                subr_loc = etree.SubElement(subcellular_loc, 'subcellularLocation')
                                sub = sc.strip().split(',')
                                for s in sub:
                                    location = etree.SubElement(subr_loc, 'location')
                                    location.text = s.strip()
                            else:
                                etex = etree.SubElement(subcellular_loc,'text')
                                etex_re = re.search('Note=(.*)', sc.strip())
                                etex.text = etex_re.group(1).strip('\n').strip('\r')
                if go_found == True:
                    for g1, g2 in zip(goid_list,goterm_list):
                        dbre_at = {'type':'GO', 'id':g1}
                        dbReference = etree.SubElement(entry, 'dbReference', dbre_at)
                        dbrep = {'type':'term', 'value':g2}
                        dbprop = etree.SubElement(dbReference, 'property', dbrep)
                if ft_found == True:
                    for ft in FT_array:
                        ft_dict = dict()
                        ft_dict = {'type':ft[0]}
                        if not ft[3] == '':
                            ft_dict['description'] = ft[3]
                        feature = etree.SubElement(entry, 'feature', ft_dict)
                        location = etree.SubElement(feature, 'location')
                        
                        if ft[1] == ft[2]:
                            posit = etree.SubElement(location, 'position', {'position': ft[1]})
                            
                        else:
                            posbe = etree.SubElement(location, 'begin', {'position': ft[1]})
                            posen = etree.SubElement(location, 'end', {'position': ft[2]})
                if seq_loc == True:
                    sequence = etree.SubElement(entry, 'sequence', {'length': seq_aa, 'mass': seq_mw})
                    #print(acc_list)
                    #print(etree.tostring(entry, pretty_print=True))
                    sequence.text = seq
                count += 1
                print(count)
                start_entry = False
                seq_loc = False
                subcc = False
                go_found = False
                ft_found = False
                acc_list = []
                subcell_loc = []
                goid_list = []
                goterm_list = []
                fn = ''
                ecn = ''
                FT_array = []
                
                
                outputfile.write(etree.tostring(entry, encoding='utf-8', pretty_print=True))
        outputfile(b'</uniprot>')   
                        
                        
               
                        
                    
                    
                    
        
                
