#Angel SNP script
from  blast_output_reader_ACI import run_RNAstructure_fold, store_rando_files, run_spot_wt_and_snp, add_snps_nomenclature, run_probknot_general, run_RNAsnp, set_up_intermediate_SNPfold, organize_genome_folder
import os 
import subprocess
import time
#from openpyxl import load_workbook
#from openpyxl import Workbook


lengths_chr = [42787722, 37316900, 38991741, 34943064, 30810063, 31706320, 29491485, 29456214, 24344104, 24739695, 31893087, 26473654]

def chr_length_to_append(chr):
    chr = int(chr)
    length = 0
    for addition in lengths_chr[:chr-1]:
        length += addition
    return length


class Mutant:
    def __init__(self, chr, position, id, wt_base, mu_base):
        self.chr = chr
        self.position = position
        self.id = id
        self.wt_base = wt_base
        self.mu_base = mu_base

    def extol(self):
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.chr, self.position, self.id,self.locus_id, self.description, self.wt_base,self.mu_base, self.effect))

class mRNA:
    def __init__(self, chr, position, id, sequence):
        self.chr = chr
        self.position = int(position)
        self.id = id
        self.subset = 'full_gene'
        self.sequence = sequence
        self.list_snps = []
        self.list_mutants = []
        self.length = len(self.sequence)
        self.end_position = int(position) + int(self.length)

class SNP:
    def __init__(self, wt, alt, n22_pos, nb_pos):
        self.wt = wt
        self.alt = alt
        self.n22_pos = n22_pos
        self.nb_pos = nb_pos
        self.snp_dg = ''
        self.wt_dg = ''
        self.ddg = ''
        self.pk_results = ''
        self.rnasnp_results = ''
        self.mrna = ''

    def initialize_list(self):
        self.list = [self.n22_pos,self.wt, self.alt, self.nb_pos, self.snp_dg, self.wt_dg, self.ddg, self.pk_results, self.rnasnp_results]

    def extol(self):
        self.initialize_list()
        new_list = []
        for i in self.list:
            new_list.append(str(i))
        extol = '\t'.join(new_list)
        return extol
#this function writes the SNP files for RNAsnp to input
def write_snp_file(wt_base, position, mu_base):
    snp_string = wt_base+str(position)+mu_base
    snp_SNP_filename = snp_string+'.fasta'
    snp_SNP_file = open(snp_SNP_filename, 'w')
    snp_SNP_file.write(snp_string)
    snp_SNP_file.close()
    return [snp_string, snp_SNP_filename]

def master():
    print(time.ctime(time.time()))
    number_runs = open('run_ID_N22.txt', 'r')
    number_runs_lines = number_runs.readlines()
    num_prev_run = int(number_runs_lines[0].replace('\n', ''))
    this_run_num = num_prev_run + 1
    number_runs.close()
    number_runs = open('run_ID_N22.txt', 'w')
    number_runs.write(str(this_run_num)+'\n')
    number_runs.close()
    #extracts all useful information from excel, stores in objects in the list
    #list_snps = information_from_excel('SNPs_from_1001G.xlsx')
    list_genes = genes_wanted('differentially_expressed_genes_n22.txt')
    list_sequences = information_from_ensemble('ensembl_gtf_ALL.fa', list_genes)
    try:
        os.mkdir('nipponbare_n22_angel')
    except:
        pass
    os.chdir('nipponbare_n22_angel')
    dir_name = 'riboSNitch_programs_data'
    try:
        os.mkdir(dir_name)
    except:
        pass
    all_snps = run_blastn_mRNA(list_sequences)
    all_files = all_snps[1]
    all_snps = all_snps[0]
    #for every SNP mutant a mutant sequence is contructed and stored
    print('\nConstructing Mutants and calculating free energy changes...\n')
    for mrna in list_sequences:
        print('\n============'+mrna.id)
        ref_filename = 'N22_Chr%s.fasta' % mrna.chr
        print(mrna.list_snps, 'KKKKKKKKKKKK')
        for entry in mrna.list_snps:
            print('\n-----------'+entry.wt+str(entry.nb_pos)+entry.alt)
            construction = construct_mutant(mrna.id+'.fasta', mrna.chr, entry.nb_pos, entry.wt, entry.alt, mrna.id+entry.wt+str(entry.nb_pos)+entry.alt)
            if construction == 'insertion' or  construction == 'deletion':
                continue
            temp_mutant = Mutant(mrna.chr, construction[2], mrna.id+entry.wt+str(entry.nb_pos)+entry.alt, entry.wt, entry.alt)
            temp_mutant.sequence = construction[0]
            temp_mutant.wt_sequence = construction[1]
            mrna.list_mutants.append(temp_mutant)
            print('Phase 1')
            for snp in all_snps:
                if entry.n22_pos == snp.n22_pos and entry.wt == snp.wt and entry.nb_pos == snp.nb_pos and entry.alt == snp.alt:
                    snp.snp_dg = (run_RNAstructure_fold(temp_mutant.id+'.fasta', snp))
                    snp.wt_dg = (run_RNAstructure_fold('ref_'+temp_mutant.id+'.fasta', snp))
                    snp_fe = float(snp.snp_dg.split(' ')[0])
                    snp_fe_er = float(snp.snp_dg.split(' ')[2])
                    wt_fe = float(snp.wt_dg.split(' ')[0])
                    wt_fe_er = float(snp.wt_dg.split(' ')[2])
                    total_err = round(wt_fe_er + snp_fe_er,2)
                    d_fe = round(snp_fe - wt_fe, 2)
                    snp.ddg = (str(d_fe)+' '+snp.snp_dg.split(' ')[1]+' '+str(total_err))
    #runs through all the SNP mutants
    for mrna in list_sequences:
        for mutant in mrna.list_mutants:
            if mutant.wt_base == mutant.mu_base:
                continue
            print(mutant.id)
            snp_filename = write_snp_file(mutant.wt_base, mutant.position, mutant.mu_base)
            #subprocess.run('RNAsnp -f %s -s %s -w 100' % ('ref_'+entry.id+ '.fasta', snp_filename[1]), shell=True)
            #runs spot on the wt and the snp sequence and compares them
            print(snp_filename, 'MMMMM', mutant.id, mutant.id, snp_filename[0], mutant.sequence, mutant.wt_sequence, mutant.id+'.fasta', 'nipponbare_n22_angel/')
            print(os.getcwd())
            comparison_spot = run_spot_wt_and_snp(mutant.id, mutant.id, snp_filename[0], mutant.sequence, mutant.wt_sequence, mutant.id+'.fasta', 'nipponbare_n22_angel/')
            track = comparison_spot[1]
            comparison = comparison_spot[0]
            
            track_pk_bp_sl_bp(comparison)
            #runs probknot on both the wt and the snp and compares them
            print(os.getcwd())
            prob_results = run_probknot_general(mutant.id, snp_filename[0], mutant.id+'.fasta','ref_'+mutant.id+'.fasta', mutant.id+'.fasta')
            print(prob_results, track, 'lllll')
            num_pk_prgs = 0
            pk_prgs = ''
            if track != None:
                num_pk_prgs += 1
                pk_prgs = track
            if prob_results != None:
                num_pk_prgs += 1
                pk_prgs = prob_results
            if track != prob_results and num_pk_prgs == 2:
                pk_prgs = 'change'
            pk_string = str(num_pk_prgs)+' '+str(pk_prgs)
            
            #runs RNAsnp
            os.chdir('../')    
            rnasnp_results = run_RNAsnp(mutant.id, snp_filename[0],'ref_'+mutant.id+'.fasta', snp_filename[1], 0, mutant.id+'.fasta', dir_name, False)
            if rnasnp_results != None:
                rnasnp_results = rnasnp_results[1]
                if rnasnp_results == 2:
                    rnasnp_results = 1

            print(pk_string, rnasnp_results, 'results')
            for place, snp in enumerate(all_snps):
                occurence = 1
                for other_place, other_snps in enumerate(all_snps):
                    if snp.wt == other_snps.wt and snp.alt == other_snps.alt and snp.n22_pos == other_snps.n22_pos and place != other_place:
                        del all_snps[other_place]
                        occurence += 1
                snp.occurence = (occurence)
            for snp in all_snps:
                print(mutant.wt_base, snp.wt, mutant.mu_base, snp.alt, str(mutant.position), str(snp.n22_pos), str(snp.nb_pos), 'test')
                print(snp)
                print(mutant.id[14:])
                snp.initialize_list()
                new_snp = snp.list
                snp_script_mini = snp.wt+str(snp.nb_pos) + snp.alt
                mu_snp_script_mini = mutant.id[14:]
                if mu_snp_script_mini == snp_script_mini:
                    print('appended')
                    snp.rnasnp_results = (rnasnp_results)
                    snp.pk_results = (pk_string)
                    set_up_intermediate_SNPfold_mod(mutant.id,snp_filename[0], 'ref_'+mutant.id+'.fasta', mutant.id+'.fasta', 0, this_run_num, new_snp, dir_name, '../../nipponbare_n22_angel/', 'False')
            #gets everything ready to run for SNPfold
    data_all_snps(all_snps, all_files, this_run_num)
    store_rando_files(all_files)
    #compile_ROI_data_mRNA(this_run_num, list_sequences)
    directory_cleaner()
    print(os.getcwd(),'line164')
    organize_genome_folder_mod()
    print(time.ctime(time.time()))

def directory_cleaner():
    for file in os.listdir(os.getcwd()):
        if '.fasta' in file and 'N22' not in file and 'ref' not in file:
            os.remove(file)

def track_pk_bp_sl_bp(comparison):
    tracking = open('pk_and_sl_bp.txt', 'a')
    tracking.write(str(comparison[1])+'\t'+ str(comparison[5])+'\n')

#this function serves to create mutant sequences for the SNPs
def construct_mutant(ref_file_name, input_chr, loci, wt_base, mu_base, mu_ID):
    #this is the reference sequence fasta file
    ref_file = open(ref_file_name, 'r')
    mutant_file_name = mu_ID+'.fasta'
    mu_file = open(mutant_file_name, 'w')
    relevant_ref_seq_name = 'ref_'+mu_ID+'.fasta'
    rel_ref_file = open(relevant_ref_seq_name, 'w')
    wt_lines = ref_file.readlines()
    input_chr = int(input_chr)
    chr_adjustment = chr_length_to_append(input_chr)
    loci = int(loci)
    #this loops through all the lines in the wt sequence
    for num_line, line in enumerate(wt_lines):
        select_chr = False
        #selects for informative lines, not sequence lines
        if line[0] == '>':
            select_chr = True
        #only allows the correct chromosome through
        if select_chr == True:
            wt_chr_seq = '' 
            for chr_line in wt_lines[num_line+1:]:
                if chr_line[0] == '>':
                    break
                wt_chr_seq += chr_line.replace('\n', '')

            mu_seq = ''
            #this ensures that the infomration about the SNP is corret before proceeding
            loci = int(loci)
            wt_base = str(wt_base)
            mu_base = str(mu_base)
            location_in_seq = 61
            location_top = 60
            if 'DEL' in wt_base:
                return 'deletion'
            if '-' in wt_base or '-' in mu_base:
                return 'insertion'
            if str(wt_base) == str(wt_chr_seq[int(loci)-1]):
                #this splices in the mutant base pair 
                mu_seq = wt_chr_seq[:loci-1] +'|' +mu_base+'|' + wt_chr_seq[loci:]
                #this selects the relevant and usable portion of the sequence that we need
                if loci >= 61:
                    location_in_seq = 61
                    mu_relevant_seq = mu_seq.replace('|','')[loci-61:loci+60]
                elif loci < 61:
                    location_in_seq = loci
                    mu_relevant_seq = mu_seq.replace('|','')[:loci+60]
                elif len(mu_seq) - 60 < loci:
                    location_top = len(mu_seq)-1
                    mu_relevant_seq = mu_seq.replace('|','')[loci-61:location_top]

                #the following stores the information
                mu_file.write('>%s\n%s'%(mutant_file_name, mu_relevant_seq))
                rel_ref_file.write('>%s\n%s'%(relevant_ref_seq_name, wt_chr_seq[loci-location_in_seq:loci+location_top]))
                #this returns relevent sequence information
                return [mu_relevant_seq, wt_chr_seq[loci-location_in_seq:loci+location_top], location_in_seq]
            else:
                print('FAIL', mutant_file_name)
   
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
    for line in text_file_lines[2:]:
        data = line.split('\t')
        if data != ['', '', '', '', '', '', '', '\n']:
            extract_mutant = Mutant(data[0], data[1],data[2],data[3],data[4],data[5], data[6],data[7])
            list_snps.append(extract_mutant)
    return list_snps
   #     col += 1
#construct_mutant('pseudo6909.fasta',3,104,'T','G',0)
def genes_wanted(genes_of_interest_filename):
    gene_file = open(genes_of_interest_filename, 'r')
    genes_lines = gene_file.readlines()
    list_genes = []
    for line in genes_lines[2:]:
        data = line.split('\t')
        if line != ['', '\n']:
            for gene in data:
                gene = gene.replace('g','t').replace('\n','').replace('\r','').replace('-', '-0')
                if gene != '':
                    print(gene, 'IIIIIII')
                    list_genes.append(gene[:-2])
    return list_genes

def information_from_ensemble(text_file_name, list_genes):
    text_file = open(text_file_name, 'r')
    gene_file_thing = open('gene_accounting.txt', 'a')
    gene_file_thing.write('Account of Already run genes:\n')
    gene_file_thing.close()
    text_file_lines = text_file.readlines()
    list_sequences = []
    gene_sequence = ''
    line_id = ''
    coding = False
    for line in text_file_lines:
        if len(list_sequences) > 50:
            print('!!!!!!')
            break
        if coding == True and line[0] != '>':
            gene_sequence += line.replace('\n','')
            print('11111')
        if line[0] == '>' and gene_sequence == '':
            coding = True
            line_id = line.replace('\n','')
            print('22222')
        if line[0] == '>' and gene_sequence != '':
            try:
                int(line_id[3:5])+1
                print('33333')
                check = checking_in_accounting_gene_file(line_id)
                print(check)
                print(line_id[1:16])
                if (check == True) and (line_id[1:16] in list_genes):
                    accounting_for_scanned_genes(line_id)
                    print(line_id, 'QQQQQ')
                    list_sequences.append(mRNA(line_id[3:5], line_id[6:13], line_id[1:].replace('\n','').replace('-',''), gene_sequence))
                elif (check == True) and (line_id[1:14] in list_genes):
                    accounting_for_scanned_genes(line_id)
                    print(line_id, 'QQQQQ2222')
                    list_sequences.append(mRNA(line_id[3:5], line_id[6:13], line_id[1:].replace('\n','').replace('-',''), gene_sequence))
            except:
                pass
            gene_sequence = ''
            line_id = line.replace('\n','') 
    return list_sequences

def checking_in_accounting_gene_file(input_line):
    print('@@@@@', input_line, os.getcwd())
    accounting = open('gene_accounting.txt', 'r')
    acc_lines = accounting.readlines()
    temp = input_line.replace('-','')+'\n'

    if temp not in acc_lines and temp != '':
        print(temp, 'OOOOOO')
        accounting.close()
        return True
    else:
        accounting.close()
        print('LLLLLL')
        return False

def accounting_for_scanned_genes(line):
    accounting = open('gene_accounting.txt', 'a')
    line_id = line.replace('\n','').replace('-','')
    if line_id[0] == '>':
        accounting.write(line_id+'\n')
        accounting.close()

def make_regular_fasta_from_mRNA(mrna):
    fasta_name = mrna.id +'.fasta'
    fasta = open(fasta_name, 'w')
    fasta.write('>'+fasta_name+'\n'+mrna.sequence)
    fasta.close()
    return fasta_name

def run_blastn_mRNA(input_mRNA_list):
    print(time.ctime(time.time()))
    all_snps = []
    all_files = []

    #this ensures that it is a genomic file that is being scanned
    #if file == 'MT451283.1Severeacuterespiratorysyndromecoronavirus.fasta':
    print('Blasting the mRNAs...')
    for mrna in input_mRNA_list:
        print('\n==========================================================\n'+mrna.id+'\n')
        ref_file = 'N22_Chr'+str(mrna.chr)+'.fasta'
        mrna_fasta = make_regular_fasta_from_mRNA(mrna)
        all_files.append(mrna.id+'.fasta')
        snp_file_file_name = mrna_fasta.replace('.fasta', '_snps.txt')
        snp_file_file = open(snp_file_file_name, 'w')
        snp_file_string = ''
        temp_file_name = mrna_fasta.replace('.fasta','_blastn.txt')
        temp_file = open(temp_file_name, 'w')
        temp_file.close()
        #this blasts the sequence
        #individual_blast(file, ref_file, temp_file_name)
        #-outfmt 1 -num_alignments 1 -num_descriptions 1
        subprocess.run('blastn -query %s -subject %s -outfmt 1 -max_hsps 1 -out %s' % (mrna_fasta, ref_file, temp_file_name),shell = True)
        #this runs the rest of the script
        file_snps = extract_data_from_blastn_mod(temp_file_name)
        mrna.list_snps = file_snps
        #this collects the sequence specific SNPs
        for lst in file_snps:
            all_snps.append(lst)
            snp_file_string += (lst.wt +' '+str(lst.nb_pos)+' '+str(lst.alt)+'\n')
        snp_file_file.write(mrna_fasta+'\n'+snp_file_string)
    return [all_snps, all_files]

def data_all_snps(all_snps, all_files, this_run_num):
    #this collects all the snps from the different sequences and makes sure there are no repeats and counts the number of specific repeats there are

    time_id = (time.ctime(time.time())).replace(' ', '_').replace(':', '.')
    all_snps_file = open('all_snps_%s_%s.txt' % (str(this_run_num), time_id), 'w')
    all_snps_string = 'wt\tloci\tmu\tOccurences\tdG SNP\tdG wt\tddG\trbSN\tpsSN\tloci_nb\tmrna\n'
    all_snps_name_string = ''
    all_snps.sort(key=lambda x:x.n22_pos)
    print('all_snps')
    for snp in all_snps:
        print(snp)
    for snp in all_snps:
        all_snps_string += ('\t'.join([snp.wt, str(snp.n22_pos), snp.alt, str(snp.occurence), str(snp.snp_dg), str(snp.wt_dg), str(snp.ddg), str(snp.rnasnp_results), str(snp.pk_results), str(snp.nb_pos), str(snp.mrna)])+'\n')
        #try: 
        #    all_snps_string += ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (snp[1], snp[0], snp[3], str(snp[10]), str(snp[5]), str(snp[6]), str(snp[7]), str(snp[8]), str(snp[9])))
        #    #all_snps_string += (snp[1]+'\t'+str(snp[0])+'\t'+snp[3]+ '\tOccurences: '+str(snp[8])+'\tdG SNP: '+str(snp[5])+'\tdG WT: '+str(snp[6])+'\tddG: '+str(snp[7])+'\n')
        #except:
        #    try:
        #        all_snps_string += ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (snp[1], snp[0], snp[3], str(snp[8]), str(snp[5]), str(snp[6]), str(snp[7]), '', ''))
        #    except:
        #        all_snps_string += ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (snp[1], snp[0], snp[3], str(snp[5]), '', '', '', '', ''))
        #    #all_snps_string += (snp[1]+'\t'+str(snp[2])+'\t'+snp[3]+ '\tOccurences: '+str(snp[5])+'\n')
    file_string = ''
    for file in all_files:
        file_string+=file+'\n'
    #for name in all_files:
        #all_snps_name_string += name+', '
    all_snps_file.write(time_id+'\n\n'+file_string+'\n\n'+all_snps_string)
    #all_snps_file.close()

def extract_data_from_blastn_mod(input_file):
    file = open(input_file, 'r')
    file_lines = file.readlines()
    list_location_snp = []
    temp_location_snp = []
    line_beginning_loci = 0
    chr_num = int(input_file.split('Os')[1][:2])
    chr_positiion = chr_length_to_append(chr_num)
    #this for loop examines each line of the data file
    first_subj = False
    for number_line, line in enumerate(file_lines):
        addition = False
        line_items = line.split(' ')
        #this conditional selects for the lines of interest
        if line_items[0] == 'Query_1':
            first_subj = True
        if line_items[0] == 'Subject_1' and first_subj == True:
            first_subj = False
            #this for loop goes thorguh each entry in the data line
            for num_item, item in enumerate(line_items):
                #this defines where in the genome the SNP is located
                #this finds and identifies the wt base

                if str(item).isdigit() == True:
                    loci_genomic = int(item)
                latter_loci_genomic = int(line_items[-1].replace('\n',''))

                if '.' in item:
                    for num_prev, prev_item in enumerate(file_lines[number_line-1].split(' ')):
                        if str(prev_item).isdigit() == True:
                            #this is the line beginning of the NB gene seqences (small number)
                            line_beginning_loci = int(prev_item)
                            break
                    for number, character in enumerate(item):
                        if character in ['A', 'G', 'C', 'T', '-']:
                                #output_genomic_loci is the final determination of where the SNP occurs in N22
                                if latter_loci_genomic < loci_genomic:
                                    output_genomic_loci = loci_genomic-number
                                else:
                                    output_genomic_loci = loci_genomic+number
                                #snp_loci is the final location of the SNP as it relates to NB == good for determining subset
                                snp_loci = number + line_beginning_loci
                                #Where the mutant is in the line of data, it's char (ALT), 
                                temp_location_snp = [number, character, snp_loci, line_beginning_loci, output_genomic_loci]
                                addition = True
                                if addition == True:
                                    prev_items = file_lines[number_line-1].split(' ')
                                    for item in prev_items:
                                        if ('A' in item) or ('T' in item) or ('C' in item) or ('G' in item):
                                            #for snpoly in temp_location_snp:
                                            snpoly = temp_location_snp
                                            for number, character in enumerate(item):
                                                if number == snpoly[0]:
                                                    ####NOTE THE switch from the original!!!
                                                    #final loci of N22, WT, line_beginning loci, ALT,
                                                    temp_snp = SNP(character.capitalize(), snpoly[1], int(snpoly[4])+chr_positiion, snpoly[2])
                                                    temp_snp.mrna = input_file.split('_')[0]
                                                    list_location_snp.append(temp_snp)
                                            break
                                #print(temp_location_snp)
                    #this loop detects insertions of the subject sequence
#                    for num_item2, item2 in enumerate(file_lines[number_line+1].split(' ')):
#                        if '\\' in item2:
#                            lstr = file_lines[number_line+3].split(' ')
#                            #' '.join(lstr).split()
#                            for num_item3, item3 in enumerate(lstr):
 #                               item3 = item3.replace('\n', '')
#                                if item3 in ['A', 'C', 'G', 'T']:
 #                                   deletion = True
#                                    position_tracker = -1
#                                    please_break  = False
#                                    for final_item in line_items:
#                                        if please_break == True:
#                                            break
#                                        if '.' not in final_item:
#                                            position_tracker += 1
#                                        for i in final_item:
#                                            if i == '.':
#                                                please_break = True
#                                                break
#                                            position_tracker += 1
#                                    location = num_item3-position_tracker
#                                    if [location, 'IN '+item3, loci_genomic+location, 'in',line_beginning_loci] not in list_location_snp:
#                                        list_location_snp.append([location, 'IN '+item3, loci_genomic+location, 'in',line_beginning_loci+location])

        #this functions to determine what the mutant base was for
        #  the SNP
    add_snps_nomenclature_mod(list_location_snp)
    return list_location_snp

def add_snps_nomenclature_mod(list_snps):
    list_names = ['R', 'Y', 'M', 'K', 'S','W','H','B','V','D','N']
    list_replacements = [['G', 'A'], ['T', 'C'], ['A', 'C'], ['G', 'T'], ['G', 'C'],['A', 'T'],['A','C', 'T'], ['G', 'T', 'C'], ['G','C', 'A'],['G','T','A'], ['G','C','A','T']]
    i = 0
    while i <= len(list_snps)-1:
        if list_snps[i].wt.capitalize() in list_names:
            index_value = list_names.index(list_snps[i].wt.capitalize())
            list_snps[i].wt = list_replacements[index_value][0]
            for replacement in list_replacements[index_value][1:]:
                temp_snp = SNP(list_snps[i].wt, replacement, list_snps[i].n22_pos,list_snps.nb_pos)
                list_snps.append(temp_snp)
        i += 1

home_directory_path = '/storage/work/p/pcf5065/'

def compile_ROI_data_mRNA(this_run_num, list_ROI, output_dir = 'nipponbare_n22_angel'):
    print('here')
    os.chdir(home_directory_path + output_dir)
    for snp_file in os.listdir(os.getcwd()):
        if 'all_snps' in snp_file:
            find_ID = snp_file.split('all_snps_')
            number_ID = find_ID[1].split('_')
            print(number_ID, this_run_num, 'HHHH')
            if number_ID[0].isdigit() == False:
                continue
            if this_run_num == int(number_ID[0]):
                print('found all_snps')
                snp_file_open = open(snp_file, 'r')
                snp_file_lines = snp_file_open.readlines()
                for ROI in list_ROI:
                    ROI_snp_lines = []
                    number_snps = 0
                    #this tracks and makes sure uncertainties don't triple the number of SNPs in an ROI
                    snp_ID = []
                    number_rbSN = 0
                    number_psSN = 0
                    for line in snp_file_lines[2:]:
                        data = line.split('\t')
                        if len(data) > 1:
                            if int(data[1]) > int(ROI.end_position):
                                break
                            elif int(ROI.position) < int(data[1]) < int(ROI.end_position):
                                ROI_snp_lines.append(line)
                                if int(data[1]) not in snp_ID:
                                    number_snps += int(data[3].split(' ')[1])
                                    snp_ID.append(int(data[1]))
                                    if len(data) > 6:
                                        if int(data[7].split(' ')[1]) != 0:
                                            number_rbSN += 1
                                        if int(data[8].split(' ')[1]) != 0:
                                            number_psSN += 1
                    print(ROI_snp_lines, 'LLLLLL')
                    ROI_dir = 'ROI_'+str(ROI.id)+'_'+str(ROI.subset)+'_Pass.drct/'
                    ROI_info_file_name = 'ROI_'+str(ROI.id)+'_'+str(ROI.subset)+'_info.txt'
                    try:
                        print('making directory')
                        os.mkdir(ROI_dir)
                        os.chdir(ROI_dir)
                        ROI_info_file = ROI_info_file = open(ROI_info_file_name, 'a')
                        ROI_snp_string = ''
                        for snp in ROI_snp_lines:
                            ROI_snp_string += snp
                        print(ROI_snp_string, 'iiiiii')
                        ROI_info_file.write(str(ROI.position)+'\t'+str(ROI.end_position) + '\n'+ROI.sequence+'\n'+ROI.id +'_'+ROI.subset+'\nNumber SNPs: '+str(number_snps)+'\nNumber rbSN: '+str(number_rbSN)+'\nNumber psSN: '+str(number_psSN)+'\n'+ROI_snp_string)
                        os.chdir(home_directory_path + output_dir)
                    except:
                        print('except')
                        os.chdir(ROI_dir)
                        ROI_info_file = open(ROI_info_file_name, 'r')
                        ROI_info_lines = ROI_info_file.readlines()
                        print(ROI_info_lines)
                        old_ROI_snp_lines = []
                        for line in ROI_info_lines[6:]:
                            old_ROI_snp_lines.append(line)

                        for line in old_ROI_snp_lines:
                            ROI_snp_lines.append(line)

                        print(ROI_snp_lines, 'KKKKKKK')
                        tracking_lines = []
                        snps_for_deletion = []
                        for place, snp in enumerate(ROI_snp_lines):
                            occurence = int(snp.split('\t')[3].split(' ')[1])
                            if snp[1:4] not in tracking_lines:
                                tracking_lines.append(snp[1:4])
                                for other_place, other_snp in enumerate(ROI_snp_lines):
                                    if snp.split('\t')[:3] == other_snp.split('\t')[:3] and place != other_place:
                                        occurence += int(other_snp.split('\t').split(' ')[1])
                                        snps_for_deletion.append(other_place)
                            data_snp = snp.split('\t')
                            ROI_snp_lines[place] = '\t'.join(data_snp[:3]) + '\t'+'Occurences: '+str(occurence)+ '\t'+'\t'.join(data_snp[4:])
                            ROI_snp_lines.sort(key=lambda x:x[2])

                        i = 1
                        while i <= len(snps_for_deletion):
                            del ROI_snp_lines[snps_for_deletion[-i]]
                            i += 1

                        number_snps = 0
                        #this tracks and makes sure uncertainties don't triple the number of SNPs in an ROI
                        snp_ID = []
                        number_rbSN = 0
                        number_psSN = 0
                        for line in ROI_snp_lines:
                            data = line.split('\t')
                            if data[1] > ROI.end_position:
                                break
                            elif ROI.position < data[1] < ROI.end_position:
                                if data[1] not in snp_ID:
                                    number_snps += int(data[3].split(' ')[1])
                                    snp_ID.append(data[1])
                                    if len(data) > 6:
                                        if int(data[7].split(' ')[1]) != 0:
                                            number_rbSN += 1
                                        if int(data[8].split(' ')[1]) != 0:
                                            number_psSN += 1
                        ROI_snp_string = ''
                        for snp in ROI_snp_lines:
                            ROI_snp_string += snp
                        print(ROI_snp_string, 'jjjjjjj')
                        ROI_info_file.close()
                        ROI_info_file = open(ROI_info_file_name, 'w')    
                        ROI_info_file.write(str(ROI.position)+'\t'+str(ROI.end_position) + '\n'+ROI.sequence+'\n'+ROI.id +'_'+ROI.subset+'\nNumber SNPs: '+str(number_snps)+'\nNumber rbSN: '+str(number_rbSN)+'\nNumber psSN: '+str(number_psSN)+'\n'+ROI_snp_string)
                        os.chdir(home_directory_path + output_dir)

def set_up_intermediate_SNPfold_mod(ID, snp_string, wt_filename, input_file, degree_support,this_run_num, snp_list, dir_name, working_dir = '', default=True):
    intermediate_snpfold_storage_file = open('intermediate.txt', 'a')
    intermediate_snpfold_storage_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n '%(str((ID)), str(snp_string), wt_filename, input_file, str(0), str(this_run_num), str(snp_list), dir_name, str(default)))
    intermediate_snpfold_storage_file.close()

def organize_genome_folder_mod():
    try:
        os.mkdir('blastn_outputs_Pass')
    except:
        pass
    try:
        os.mkdir('collected_snps_Pass')
    except:
        pass

    for file in os.listdir():
        if 'Pass' not in file:
            if 'blastn' in file:
                os.rename(file, 'blastn_outputs_Pass/'+file)
            if 'snp' in file:
                os.rename(file, 'collected_snps_Pass/'+file)

master()