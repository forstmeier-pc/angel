
from angel_snp_arabdopsis import construct_mutant
import subprocess
import os

lengths_chr = [42787722, 37316900, 38991741, 34943064, 30810063, 31706320, 29491485, 29456214, 24344104, 24739695, 31893087, 26473654]

def chr_length_to_append(chr):
    chr = int(chr)
    length = 0
    for addition in lengths_chr[:chr-1]:
        length += addition
    return length

def get_snps(snp_file):
    s = open(snp_file,'r')
    output_lines_list = []
    lines = s.readlines()
    output_lines_list= [lines[0].replace('\n', ',RiboSNitch\n')]
    for line in lines[1:]:
        output_lines_list.append(line)
        riboSNitch = False
        data = line.split(',')
        if len(data) < 2:
            continue
        chromo = data[0].replace('\n','')
        positi = data[1].replace('\n','')
        wt_ref = data[6].replace('\n','')
        alt_ref = data[7].replace('\n','')
        ID = data[0]+'_'+str(positi)+'_'+str(chromo)
        if len(wt_ref) > 1 or len(alt_ref) > 1:
            continue
        op = construct_mutant('TAIR10_chr_all.fas',chromo,positi,wt_ref,alt_ref,ID)
        
        if op != None:
            SNPfold_output_name = run_SNPfold_2(ID, wt_ref, positi,alt_ref,op[0])
            if SNPfold_output_name == None:
                continue
            subprocess.run('head '+SNPfold_output_name, shell=True)
            sf = open(SNPfold_output_name,'r')
            sf_data = sf.readlines()[1].split('\t')
            if float(sf_data[1]) <= 0.8 and float(sf_data[2].replace('\n','')) <= 0.8:
                output_lines_list[-1] = output_lines_list[-1].replace('\n',',True\n')
            else:
                output_lines_list[-1] = output_lines_list[-1].replace('\n',',False\n')
            os.remove(SNPfold_output_name)
        else:
            print('Nonetype FAIL '+ ID)
        
    s = open(snp_file.replace('.csv','_rbSN.csv'),'w')
    output_string = ''
    for line in output_lines_list:
        output_string += line
    s.write(output_string)



def run_SNPfold_2(input_ID, wt, pos, alt, wt_filename):
    if wt == alt:
        return None
    os.chdir('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/')
    wt_filename = '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work/'+wt_filename
    seq = open(wt_filename,'r').readlines()[1]
    #print(os.getcwd(),'64')
    new_file= open('temp%s.txt'%(input_ID+pos), 'w')
    new_file.write(seq)
    new_file.close()
    addends = write_snp_file(wt, pos, alt)
    SNPfold_output_name = 'SNPfold_output%s_%s.txt' % (addends[0], input_ID)
    SNPfold_output = open(SNPfold_output_name, 'w')

    python_2_command = 'python2.7 SNPfold_commandline.py %s %s > %s' % ('temp%s.txt'%(input_ID+pos), addends[0], SNPfold_output_name)
    subprocess.run(python_2_command, shell=True)
    os.remove('temp%s.txt'%(input_ID+pos))
    os.chdir('/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work')
    return '../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/'+SNPfold_output_name

def run_RNAsnp_2(input_ID, wt, pos, alt, wt_filename):
    if wt == alt:
        return None
    addends = write_snp_file(wt, pos, alt)
    RNAsnp_output_name = 'RNAsnp_output%s_%s.txt' % (addends[0], input_ID)
    RNAsnp_output = open(RNAsnp_output_name, 'w')
    #print('RNAsnp -f %s -s %s -w 100 > %s' % (wt_filename, addends[1], RNAsnp_output_name))
    subprocess.run('RNAsnp -f %s -s %s -w 100 > %s' % (wt_filename, addends[1], RNAsnp_output_name), shell=True)
    return RNAsnp_output_name

def write_snp_file(wt_base, position, mu_base):
    snp_string = wt_base+str(61)+mu_base
    snp_SNP_filename = snp_string+'.snp'
    snp_SNP_file = open(snp_SNP_filename, 'w')
    snp_SNP_file.write(snp_string)
    snp_SNP_file.close()
    return [snp_string, snp_SNP_filename]

def construct_mutant(ref_file_name, input_chr, loci, wt_base, mu_base, mu_ID):
    #this is the reference sequence fasta file

    ref_file = open(ref_file_name, 'r')
    mutant_file_name = mu_ID+'.fasta'
    mu_file = open(mutant_file_name, 'w')
    relevant_ref_seq_name = 'ref_'+mu_ID+'.fasta'
    rel_ref_file = open(relevant_ref_seq_name, 'w')
    wt_lines = ref_file.readlines()
    chrom_lines = []
    #this loops through all the lines in the wt sequence
    for num_line, line in enumerate(wt_lines):
        select_chr = False
        #selects for informative lines, not sequence lines
        if line[0] == '>':
            items = line.split(' ')
            items[4] = line[1]
            if items[4].isdigit()==False or items[4] == '' or input_chr == '':
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
            mu_base = str(mu_base).replace('\n','')
            wt_base = str(wt_base).replace('\n','')
            if str(wt_base)== str(wt_chr_seq[int(loci)-1]).upper():
                #this splices in the mutant base pair 
                mu_seq = wt_chr_seq[:loci-1] +'|' +mu_base+'|' + wt_chr_seq[loci:]
                #this selects the relevant and usable portion of the sequence that we need
                mu_relevant_seq = mu_seq.replace('|','')[loci-61:loci+60]
                #the following stores the information
                mu_file.write('>%s\n%s'%(mutant_file_name, mu_relevant_seq))
                rel_ref_file.write('>%s\n%s'%(relevant_ref_seq_name, wt_chr_seq[loci-61:loci+60]))
                #this returns relevent sequence information
                ref_file.close()

                return [relevant_ref_seq_name, mutant_file_name]
            else:
                print('FAIL', mutant_file_name)

    ref_file.close()


#get_snps('1000ish_SNPs.csv')