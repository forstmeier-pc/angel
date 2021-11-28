#Angel SNP script
#from  blast_output_reader_ACI import run_spot_wt_and_snp, run_probknot_general, run_RNAsnp, set_up_intermediate_SNPfold, run_RNAstructure_fold
import os 
import subprocess
#from openpyxl import load_workbook
#from openpyxl import Workbook

class mutant:
    def __init__(self, chr, position, id, locus_id, description, wt_base, mu_base, effect):
        self.chr = chr
        self.position = position
        self.id = id.replace(' ','')
        self.locus_id = locus_id
        self.description = description
        self.wt_base = wt_base
        self.mu_base = mu_base
        self.effect = effect

    def extol(self):
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.chr, self.position, self.id,self.locus_id, self.description, self.wt_base,self.mu_base, self.effect))

#this function writes the SNP files for RNAsnp to input
def write_snp_file(wt_base, position, mu_base):
    snp_string = wt_base+str(position)+mu_base
    snp_SNP_filename = snp_string+'.fasta'
    snp_SNP_file = open(snp_SNP_filename, 'w')
    snp_SNP_file.write(snp_string)
    snp_SNP_file.close()
    return [snp_string, snp_SNP_filename]

def master():
    input_file = 'extra_SNPs_12_16_20.csv'
    #extracts all useful information from excel, stores in objects in the list
    #list_snps = information_from_excel('SNPs_from_1001G.xlsx')
    print('started')
    list_snps = information_from_text(input_file)
    print('snps gathered')
    try:
        os.mkdir('Angel')
    except:
        pass
    os.chdir('Angel')
    dir_name = 'riboSNitch_programs_data'
    try:
        os.mkdir(dir_name)
    except:
        pass
    ref_filename = 'TAIR_ref_seq_full.fasta'
    #for every SNP mutant a mutant sequence is contructed and stored
    for entry in list_snps:
        construction = construct_mutant(ref_filename, entry.chr, entry.position, entry.wt_base, entry.mu_base, entry.id)
        entry.sequence = construction[0]
        entry.wt_sequence = construction[1]
    print('contruciton complete')
    
    #small block that removes invalid SNPs from faulty data
    i = 0
    while i < len(list_snps):
        entry = list_snps[i]
        test_file = open(entry.id+'.fasta', 'a')
        test_file.write('A')
        test_file.close()
        test_file = open(entry.id+'.fasta', 'r')
        test = test_file.readlines()
        if test[0] == 'A':
            list_snps.remove(entry)
        else:
            i +=1

    #runs through all the SNP mutnats
    for entry in list_snps:
        
        print(entry.id, os.getcwd())
        snp_filename = write_snp_file(entry.wt_base, 61, entry.mu_base)
        #subprocess.run('RNAsnp -f %s -s %s -w 100' % ('ref_'+entry.id+ '.fasta', snp_filename[1]), shell=True)

        rnafold_results = both_rna_structures(entry)

        #runs spot on the wt and the snp sequence and compares them
        spot_results = run_spot_wt_and_snp(entry.id, entry.id, snp_filename[0], entry.sequence, entry.wt_sequence, entry.id+'.fasta', 'Angel/')

        #runs probknot on both the wt and the snp and compares them
        print(os.getcwd())
        print(os.getcwd())

        prob_results = run_probknot_general(entry.id, snp_filename[0], entry.id+'.fasta','ref_'+entry.id+'.fasta', entry.id+'.fasta')

        print(os.getcwd())
        
        #runs RNAsnp
        os.chdir('../')    
        print(os.getcwd())

        rnasnp_results = run_RNAsnp(entry.id, snp_filename[0],'ref_'+entry.id+'.fasta', snp_filename[1], 0, entry.id+'.fasta', dir_name, False)
        print(os.getcwd())


        #gets everything ready to run for SNPfold
        set_up_intermediate_SNPfold(entry.id,snp_filename[0], 'ref_'+entry.id+'.fasta', entry.id+'.fasta', 0, input_file,int(entry.position), dir_name, '../../Angel/', 'False')

        add_results_to_original_file(entry, spot_results, prob_results, rnasnp_results, rnafold_results, input_file)

def add_results_to_original_file(entry, spot_results, prob_results, rnasnp_results, rnafold,input_file):
    #print(str(entry.position),spot_results, prob_results, rnasnp_results, str(rnafold), 'RESULTS')
    file= open(input_file.replace('.csv','_pipeline.csv'),'r')
    lines = file.readlines()
    #print(lines, 'lines')
    if 'Prob' not in lines[0]:
        header = lines[0].replace('\n','\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%('SNP_dG','WT_dG', 'ddG','SPOT-RNA','ProbKnot','RNAsnp', 'SNPfold')).replace('\t',',')
    else:
        header = ''
    outputstring = header
    outputstring2 = header
    if spot_results != None:
        print(spot_results,'spot_results')
        spot_results = spot_results[1]
    else:
        spot_results = 'None'
    if prob_results != None:
        print(prob_results, 'prob_results')
        prob_results = prob_results
    else:
        prob_results = 'None'
    if rnasnp_results != None:
        rnasnp_results = rnasnp_results[0]
        if rnasnp_results == True:
            rnasnp_results = '1'
        elif rnasnp_results == False:
            rnasnp_results ='0'
    else:
        rnasnp_results = ''

    new_line_end = '%s,%s,%s,%s,%s,%s,%s'%(rnafold[0],rnafold[1], rnafold[2],spot_results,prob_results,rnasnp_results, 'SNPfold')
    #print(new_line_end)
    for line in lines[1:]:
        #print(line, 'line')
        output_line = line.replace('\n','')
        output_line2 = line.replace('\n','')
        data = line.split(',')
        if int(data[1]) == int(entry.position):
            print('Added new line end')
            if output_line[-1] == ',':
                output_line += new_line_end
                output_line2 += new_line_end
            else:
                output_line += (","+new_line_end)
                output_line2 += (","+new_line_end)
            outputstring2 += (output_line2+'\n')
        
        outputstring += (output_line+'\n')
    print(outputstring,'11111111')
    print(outputstring2,'2222222222')
    outputstring2.replace('\xc2','')
    file.close()
    file2 = open(input_file.replace('.csv','_pipeline.csv'),'w')
    file2.write(outputstring)
    file2.close()
    file3 = open(input_file.replace('.csv','_pipeline2.csv'),'a')
    file3.write(outputstring2)
    file3.close()

def both_rna_structures(entry):
    print('start fold')
    wt_name = 'ref_'+entry.id+'.fasta'
    wt_fasta = open(wt_name,'w')
    wt_fasta.write('>'+wt_name+'\n'+entry.wt_sequence+'\n')
    wt_fasta.close()
    print('found wt')
    snp_name = entry.id +'.fasta'
    snp_fasta = open(snp_name,'w')
    snp_fasta.write('>'+snp_name+'\n'+entry.sequence+'\n')
    snp_fasta.close()
    print('found snp')
    wt = run_RNAstructure_fold(wt_name,'')
    print('wt done')
    snp = run_RNAstructure_fold(snp_name,'')

    wt_fe = float(wt.split(' ')[0])
    wt_er = float(wt.split(' ')[2])
    snp_fe = float(snp.split(' ')[0])
    snp_er = float(snp.split(' ')[2])
    ddg = round(snp_fe-wt_fe,2)
    total_error = round(snp_er+wt_er,2)
    final_ddg_format = str(ddg)+' +/- '+str(total_error)
    snp_format = str(snp_fe)+' +/- '+str(snp_er)
    wt_format = str(wt_fe)+' +/- '+str(wt_er )
    return(snp_format,wt_format, final_ddg_format)

#this function serves to create mutant sequences for the SNPs
def construct_mutant(ref_file_name, input_chr, loci, wt_base, mu_base, mu_ID):
    #this is the reference sequence fasta file
    ref_file = open(ref_file_name, 'r')
    mutant_file_name = mu_ID+'.fasta'
    mu_file = open(mutant_file_name, 'w')
    relevant_ref_seq_name = 'ref_'+mu_ID+'.fasta'
    rel_ref_file = open(relevant_ref_seq_name, 'w')
    wt_lines = ref_file.readlines()
    #this loops through all the lines in the wt sequence
    for num_line, line in enumerate(wt_lines):
        select_chr = False
        #selects for informative lines, not sequence lines
        if line[0] == '>':
            items = line.split(' ')
            if items[4] == '' or input_chr == '':
                continue
            if int(items[4]) == int(input_chr):
                select_chr = True
        #only allows the correct chromosome through
        #print('103')
        if select_chr == True:
            wt_chr_seq = '' 
            for chr_line in wt_lines[num_line+1:]:
                if chr_line[0] == '>':
                    break
                wt_chr_seq += chr_line.replace('\n', '')

            mu_seq = ''
            #this ensures that the infomration about the SNP is corret before proceeding
            #print(str(len(wt_chr_seq)), loci)
            loci = int(loci)
            wt_base = str(wt_base)
            mu_base = str(mu_base)
            if str(wt_base) == str(wt_chr_seq[int(loci)-1]):
                #this splices in the mutant base pair 
                mu_seq = wt_chr_seq[:loci-1] +'|' +mu_base+'|' + wt_chr_seq[loci:]
                #this selects the relevant and usable portion of the sequence that we need
                mu_relevant_seq = mu_seq.replace('|','')[loci-61:loci+60]
                #the following stores the information
                mu_file.write('>%s\n%s'%(mutant_file_name, mu_relevant_seq))
                rel_ref_file.write('>%s\n%s'%(relevant_ref_seq_name, wt_chr_seq[loci-61:loci+60]))
                #this returns relevent sequence information
                ref_file.close()
                return [mu_relevant_seq, wt_chr_seq[loci-61:loci+60]]
            else:
                print('FAIL', mutant_file_name)
            #else:
            #    print(mutant_file_name, 'missed bp')
            #print (wt_chr_seq[loci-1])
            #print(wt_chr_seq[:201])
            #print(mu_seq[:201])
            #print(mu_seq.replace('|','')[:201])

    #wt_seq = wt_lines[1:]
    #print(wt_seq[2])
    #wt_seq = wt_lines[1:].replace('\n', '')
    #print(wt_seq)

#this function serves to extract the information for the excel book that Angel provided us
def information_from_excel(excel_file_name):
    excel_book = load_workbook(excel_file_name)
    sheet = excel_book.active
    colnum = 1
    rownum = 1
    list_snps = []

    while sheet.cell(row=rownum, column=1).value != None or sheet.cell(row=rownum+1, column=1).value != None:
        if sheet.cell(row=rownum, column=1).value == None or sheet.cell(row=rownum, column=1).value == 'Chr':
            rownum += 1
            continue
        #print(sheet.cell(row=rownum, column=1).value)

        #this creates objects with all the useful information in one place
        extracted_mutant = mutant(sheet.cell(row=rownum, column=1).value,sheet.cell(row=rownum, column=2).value,sheet.cell(row=rownum, column=3).value,sheet.cell(row=rownum, column=4).value,sheet.cell(row=rownum, column=5).value,sheet.cell(row=rownum, column=6).value,sheet.cell(row=rownum, column=7).value,sheet.cell(row=rownum, column=8).value)
        #extracted_mutant.extol()
        list_snps.append(extracted_mutant)
        rownum +=1
    #while col <= 8:
    #this outputs a list of the objects
    return list_snps

def information_from_text(text_file_name):
    text_file = open(text_file_name, 'r')
    text_file_lines = text_file.readlines()
    list_snps = []
    for line in text_file_lines[1:]:
        data = line.split(',')
        print(str(data),'data')
        if data != ['', '', '', '', '', '', '', '\n']:
            extract_mutant = mutant(data[0], data[1],data[2],data[3],data[4],data[5], data[6],data[7])
            list_snps.append(extract_mutant)
    print(list_snps)
    pipeline_out = ''
    for line in text_file_lines:
        pipeline_out += line
    print(pipeline_out, 'pipeline_out')
    pipeline_file = open('Angel/'+text_file_name.replace('.csv','_pipeline.csv'),'w')
    pipeline_file.write(pipeline_out)
    pipeline_file.close()
    text_file.close()
    return list_snps
   #     col += 1
#construct_mutant('pseudo6909.fasta',3,104,'T','G',0)
#master()

#construct_mutant('T')