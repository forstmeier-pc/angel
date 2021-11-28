
import subprocess
import os
import argparse

lengths_chr = [42787722, 37316900, 38991741, 34943064, 30810063, 31706320, 29491485, 29456214, 24344104, 24739695, 31893087, 26473654]

parser = argparse.ArgumentParser()

parser.add_argument('-i', type=str, help='input file')

args = parser.parse_args()

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
        #            op = construct_mutant('arabidopsis_ref_seq/GCF_000001735.4_TAIR10.1_genomic.fna',data[1],data[2],data[3],data[6],data[0]+'_'+str(data[2])+'_'+data[5].replace('\n',''))
        #if len(data) > 2:

        #    op=construct_mutant('arabidopsis_ref_seq/GCF_000001735.4_TAIR10.1_genomic.fna', data[0], data[1], data[4], data[5],str(data[2])+'_'+data[5].replace('\n',''))
        '''
            RNAsnp_output_name = run_RNAsnp_2(ID, data[6].replace('\n',''), data[2].replace('\n',''),data[3].replace('\n',''),op[0])
            if RNAsnp_output_name == None:
                continue
            r= open(RNAsnp_output_name, 'r')
            output_data = r.readlines()[1].split('\t')
            print(output_data, '1')
            print(output_lines_list[-1], '2')
            print(output_lines_list[-1].replace('\n',',True\n'),'3')
            if float(output_data[6]) <= 0.1 and float(output_data[9].replace('\n', '')) <= 0.1:
                riboSNitch = True
                output_lines_list[-1] = output_lines_list[-1].replace('\n',',True\n')
            else:
                output_lines_list[-1] = output_lines_list[-1].replace('\n',',False\n')
        '''
        
    s = open(snp_file.replace('.csv','_rbSN.csv'),'w')
    output_string = ''
    for line in output_lines_list:
        output_string += line
    s.write(output_string)

def run_SNPfold_2(input_ID, wt, pos, alt, wt_filename):
    print('82')
    if wt == alt:
        return None
    #os.chdir('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/')
    os.chdir('../programs/SNPfold-1.01/')
    print('87')
    wt_filename = '../../angel_9_14_21/'+wt_filename
    #wt_filename = '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work/'+wt_filename

    seq = open(wt_filename,'r').readlines()[1]
    print('90')
    print(input_ID)
    print(seq)
    new_file= open('temp%s.txt'%(input_ID), 'w')
    new_file.write(seq)
    new_file.close()
    print('96')
    addends = write_snp_file(wt, pos, alt)
    print('97')
    SNPfold_output_name = 'SNPfold_output%s_%s.txt' % (addends[0], input_ID)
    print('91')
    SNPfold_output = open(SNPfold_output_name, 'w')

    python_2_command = 'python2.7 SNPfold_commandline.py %s %s > %s' % ('temp%s.txt'%(input_ID), addends[0], SNPfold_output_name)
    os.system(python_2_command)
    print('101')
    os.remove('temp%s.txt'%(input_ID))
    print('102')
    os.remove(wt_filename)
    print('103')
    os.remove(wt_filename.replace('ref_',''))
    print('104')
    os.chdir('../../angel_9_14_21')
    #os.chdir('/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work')
    return '../programs/SNPfold-1.01/'+SNPfold_output_name
    #return ('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/'+SNPfold_output_name)

def run_SNPfold_2_local(input_ID, wt, pos, alt, wt_filename, length):
    if wt == alt:
        return None
    os.chdir('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/')

    wt_filename = '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work/'+wt_filename

    seq = open(wt_filename,'r').readlines()[1]
    new_file= open('temp%s.txt'%(input_ID+str(length)), 'w')
    new_file.write(seq)
    new_file.close()
    addends = write_snp_file(wt, pos, alt)
    SNPfold_output_name = 'SNPfold_output%s_%s_%s.txt' % (addends[0], input_ID, str(length))
    SNPfold_output = open(SNPfold_output_name, 'w')

    python_2_command = 'python2.7 SNPfold_commandline.py %s %s > %s' % ('temp%s.txt'%(input_ID+str(length)), addends[0], SNPfold_output_name)
    p = subprocess.run(python_2_command,shell=True)
    #P.kill()
    #print('101')
    try:
        os.remove('temp%s.txt'%(input_ID+str(length)))
    except:
        pass
    #print('102')
    os.remove(wt_filename)
    #print('103')
    os.remove(wt_filename.replace('ref_',''))
    #print('104')
    #os.chdir('../../angel_9_14_21')
    os.chdir('/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work')
    #return '../programs/SNPfold-1.01/'+SNPfold_output_name
    return ('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/'+SNPfold_output_name)

def run_SNPfold_2_local_length(input_ID, wt, pos, alt, wt_filename, length):
    print('82')
    if wt == alt:
        return None
    print('85')
    #os.chdir('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/')
    os.chdir('../programs/SNPfold-1.01/')
    print('87')
    wt_filename = '../../angel_9_14_21/'+wt_filename
    #wt_filename = '/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work/'+wt_filename
    print(wt_filename)

    seq = open(wt_filename,'r').readlines()[1]
    print('90')
    print(input_ID)
    print(seq)
    new_file= open('temp%s.txt'%(input_ID+str(length)), 'w')
    new_file.write(seq)
    new_file.close()
    print('96')
    addends = write_snp_file(wt, pos, alt)
    print('97')
    SNPfold_output_name = 'SNPfold_output%s_%s_%s.txt' % (addends[0], input_ID, str(length))
    print('91')
    SNPfold_output = open(SNPfold_output_name, 'w')

    python_2_command = 'python2.7 SNPfold_commandline.py %s %s > %s' % ('temp%s.txt'%(input_ID+str(length)), addends[0], SNPfold_output_name)
    p = os.system(python_2_command)
    #P.kill()
    print('101')
    try:
        os.remove('temp%s.txt'%(input_ID+str(length)))
    except:
        pass
    print('102')
    os.remove(wt_filename)
    print('103')
    os.remove(wt_filename.replace('ref_',''))
    print('104')
    os.chdir('../../angel_9_14_21')
    #os.chdir('/mnt/c/Users/pfors/Desktop/Research Lab/Bevilacqua Research Lab/work')
    return '../programs/SNPfold-1.01/'+SNPfold_output_name
    #return ('../../../../../../../../home/pfors/programs/SNPfold-1.01/build/scripts-3.6/'+SNPfold_output_name)

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
    snp_string = wt_base+str(int(position))+mu_base
    snp_SNP_filename = snp_string+'.snp'
    snp_SNP_file = open(snp_SNP_filename, 'w')
    snp_SNP_file.write(snp_string)
    snp_SNP_file.close()
    return [snp_string, snp_SNP_filename]

def construct_mutant(ref_file_name, input_chr, loci, wt_base, mu_base, mu_ID, constraint = 0):
    #this is the reference sequence fasta file
    ref_file = open(ref_file_name, 'r')
    mutant_file_name = mu_ID+'.fasta'
    mu_file = open(mutant_file_name, 'w')
    relevant_ref_seq_name = 'ref_'+mu_ID+'.fasta'
    rel_ref_file = open(relevant_ref_seq_name, 'w')
    wt_lines = ref_file.readlines()
    chrom_lines = []
    #this block helps construct SNPs near the terminus of the transcript
    prime5contraint = False
    prime3constraint = False
    equal = False
    if constraint != 0:
        if constraint < int(loci):
            prime5contraint = True
        if constraint > int(loci):
            prime3constraint = True
        if constraint == int(loci):
            equal = True
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
                left = loci-41
                right = loci+40
                minus = 41
                if constraint != 0:
                    if prime5contraint == True:
                        left = constraint
                        minus = loci-constraint
                    if prime3constraint == True:
                        right = constraint
                    if equal == True:
                        left = constraint-1
                        minus = 1
                #this splices in the mutant base pair 
                mu_seq = wt_chr_seq[:loci-1] +'|' +mu_base+'|' + wt_chr_seq[loci:]
                #this selects the relevant and usable portion of the sequence that we need
                mu_relevant_seq = mu_seq.replace('|','')[left:right]
                #the following stores the information
                mu_file.write('>%s\n%s'%(mutant_file_name, mu_relevant_seq))
                rel_ref_file.write('>%s\n%s'%(relevant_ref_seq_name, wt_chr_seq[left:right]))
                #this returns relevent sequence information
                ref_file.close()
                print('---')
                print(constraint)
                print(left,right,minus)
                print(wt_chr_seq[left:right])
                print(len(wt_chr_seq[left:right]))
                print('--')
                return [relevant_ref_seq_name, mutant_file_name, minus]
            else:
                print('FAIL', mutant_file_name)

    ref_file.close()

def open_gff():
    t = open('TAIR10_GFF3_genes.gff', 'r')
    lines = t.readlines()
    outlines = []
    for line in lines:
        if 'mRNA' in line or 'UTR' in line:
            outlines.append(line)
    t.close()
    return outlines

def is_SNP_near_terminus(gff_lines, chr, loci):
    for line in gff_lines:
        data = line.split('\t')
        line_chr = data[0]
        line_subset = data[2]
        line_start = int(data[3])
        line_end = int(data[4])
        if 'Chr%s' % (str(chr)) == line_chr:
            if (line_subset == 'mRNA') or ('UTR' in line_subset):
                if (within_40(int(loci), line_start) == True):
                    return line_start
                if (within_40(int(loci), line_end) == True):
                    return line_end
    return 0
            
def within_40(num, num2):
    for i in range(num-40, num+40):
        if i == num2:
            return True
    return False
#get_snps('1000ish_SNPs.csv')