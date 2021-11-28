import multiprocessing as mp
import arabidopsis_scratch as arab
import subprocess
import os
import argparse


#print("number of cpu: ", mp.cpu_count()-2)
def construct_mutant_custom(ref_file_name, input_chr, loci, wt_base, mu_base, mu_ID, upstream, downstream, constraint = 0):
    #this is the reference sequence fasta file
    print('Here now birches')
    length_length = int(upstream)+int(downstream)
    ref_file = open(ref_file_name, 'r')
    mutant_file_name = mu_ID+str(length_length)+'.fasta'
    mu_file = open(mutant_file_name, 'w')
    relevant_ref_seq_name = 'ref_'+mu_ID+str(length_length)+'.fasta'
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
                left = loci-1-downstream
                right = loci+upstream
                minus = downstream+1
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
                print(upstream, downstream)
                print(wt_chr_seq[left:right])
                print(len(wt_chr_seq[left:right]))
                #print('--')
                return [relevant_ref_seq_name, mutant_file_name, minus, len(wt_chr_seq[left:right])]
            else:
                print('FAIL', mutant_file_name)

    ref_file.close()

def work(line_data):
    #try:
    print('work')
    output_line = line_data
    data=line_data.split(',')
    pos = data[1]
    if pos.isdigit() == False or len(data)<10:
        return
    print('Process %s is working on line %s' % (mp.current_process().name, str(pos)))
    print(data)
    chromo = data[0].replace('\n','').replace('"','')
    positi = data[1].replace('\n','').replace('"','')
    info  = data[13].replace('\n','').replace('"','')
    wt_ref = data[9].replace('\n','').replace('"','')
    alt_ref = data[10].replace('\n','').replace('"','')
    i = data[len(data)-2].replace('\n','').replace('"','')
    ii = data[len(data)-1].replace('\n','').replace('"','')
    #if is_climate(chromo, positi) == 'non_climate':
    #    return
    ID = data[0]+'_'+str(positi)+'_'+str(chromo)
    print(len(wt_ref))
    if len(wt_ref) > 1 or len(alt_ref) > 1:
        return
    print(data)
    op = construct_mutant_custom('TAIR10_chr_all.fas',chromo,positi,wt_ref,alt_ref,ID,int(i),int(i)+int(ii), is_SNP_near_terminus_alt(gff_lines,chromo,positi, info))
    
    print(ID+alt_ref, str(op[3]))
    #if (ID+alt_ref)== '5_15841744_5A':
    #    return
    if op != None:
        print('26')
        SNPfold_output_name = arab.run_SNPfold_2_local_length(ID, wt_ref,op[2],alt_ref,op[0], str(op[3]))
        if SNPfold_output_name == None:
            print(SNPfold_output_name, '27')
            return
        print(ID+alt_ref)
        #os.system('head '+SNPfold_output_name)
        sf = open(SNPfold_output_name,'r')
        sf_data = sf.readlines()[1].split('\t')
        print(str(sf_data))
        sf.close()
        output_line = output_line.replace('\n',',%s,%s\n'%(sf_data[1], sf_data[2]))
        os.remove(SNPfold_output_name)
    else:
        print('Nonetype FAIL '+ ID)
    return output_line
    #except:
    #    print('Fatal Error line: '+line_data)
    

def is_SNP_near_terminus_alt(gff_lines, chr, loci, info):
    if 'UTR' not in info:
        return 0
    print(chr)
    for line in gff_lines[int(chr)]:
        data = line.split('\t')
        if len(data) < 5:
            continue
        line_chr = data[0]
        line_subset = data[2]
        line_start = int(data[3])
        line_end = int(data[4])
        if (line_subset == 'mRNA') or ('UTR' in line_subset):
            if (arab.within_40(int(loci), line_start) == True):
                return line_start
            if (arab.within_40(int(loci), line_end) == True):
                return line_end
    return 0

def pool_manager(data_lines, filename):
    #p=mp.Pool(6)
    #print('pool')
    print(mp.cpu_count())
    p = mp.Pool(mp.cpu_count())
    print('pool object')
    print(p)
    new_lines = p.map(work, data_lines)
    p.terminate()
    output_string = ''
    #print(new_lines)
    for line in new_lines:
        if line != None:
            output_string += line
    #print(output_string)
    output_filename = filename.replace('.csv','_mp_rbSN_1.csv')
    output_file = open(output_filename,'a')
    output_file.write(output_string)
    output_file.close()

def open_climate(climatename):
    c = open(climatename, 'r')
    cl = c.readlines()
    c.close()
    return cl

def is_climate(chr, pos):
    for line in climate_lines:
        data = line.split(',')
        #print(line)
        #print(data)
        if str(data[0]) == str(chr) and str(data[1]) == str(pos):
            return data[2]

def open_gffs():
    all_lines = []
    for i in range(1,6):
        t = open('Arabidopsis_thaliana.TAIR10.51.chromosome.%s.gff3'%(str(i)), 'r')
        lines = t.readlines()
        outlines = []
        for line in lines:
            if 'mRNA' in line or 'UTR' in line:
                outlines.append(line)
        t.close()
        all_lines.append(outlines)
    print(len(all_lines))
    print(len(all_lines[0]))
    return outlines


if __name__ == '__main__':
    print('ye')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input file')

    args = parser.parse_args()
    filename = args.i
    #climatename = args.c
    s = open(filename,'r')
    output_lines_list = []
    lines = s.readlines()
    data_lines = []
    gff_lines = open_gffs()
    #climate_lines = open_climate(climatename)

    for line in lines:
        data = line.split(',')
        if len(data) > 3 and data[0] != '':
            #print('---')
            #print(line)
            #print(line.replace('\n',','+str(i)+','+ str(ii)+'\n'))
            data_lines.append(line.replace('\n',','+str(40)+','+ str(0)+'\n'))
            if (len(data_lines) == mp.cpu_count()-1):
                print(len(data_lines))
                pool_manager(data_lines,filename)
                data_lines = []
