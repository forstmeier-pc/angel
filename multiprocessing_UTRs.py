import multiprocessing as mp
import arabidopsis_scratch as arab
import subprocess
import os
import argparse

#print("number of cpu: ", mp.cpu_count()-2)

def work(line_data):
    #try:
    print('work')
    output_line = ''
    data=line_data.split(',')
    pos = data[1]
    print('Process %s is working on line %s' % (mp.current_process().name, str(pos)))
    chromo = data[0].replace('\n','').replace('"','')
    positi = data[1].replace('\n','').replace('"','')
    wt_ref = data[9].replace('\n','').replace('"','')
    alt_ref = data[10].replace('\n','').replace('"','')
    ID = data[0]+'_'+str(positi)+'_'+str(chromo)
    if len(wt_ref) > 1 or len(alt_ref) > 1:
        return
    op = arab.construct_mutant('TAIR10_chr_all.fas',chromo,positi,wt_ref,alt_ref,ID, arab.is_SNP_near_terminus(gff_lines,chromo,positi))
    if op != None:
        print('26')
        SNPfold_output_name = arab.run_SNPfold_2(ID, wt_ref,op[2],alt_ref,op[0])
        if SNPfold_output_name == None:
            print(SNPfold_output_name, '27')
            return
        print(ID+alt_ref)
        #os.system('head '+SNPfold_output_name)
        sf = open(SNPfold_output_name,'r')
        sf_data = sf.readlines()[1].split('\t')
        print(str(sf_data))
        sf.close()
        if float(sf_data[1]) <= 0.8 and float(sf_data[2].replace('\n','')) <= 0.8:
            output_line = line_data.replace('\n',',True\n')
        else:
            output_line = line_data.replace('\n',',False\n')
        os.remove(SNPfold_output_name)
    else:
        print('Nonetype FAIL '+ ID)
    return output_line
    #except:
    #    print('Fatal Error line: '+line_data)


def pool_manager(data_lines, filename):
    #p=mp.Pool(1)
    print('pool')
    print(mp.cpu_count())
    p = mp.Pool(mp.cpu_count())
    print('pool object')
    print(p)
    new_lines = p.map(work, data_lines)
    output_string = ''
    print(new_lines)
    for line in new_lines:
        if line != None:
            output_string += line
    print(output_string)
    output_filename = filename.replace('.csv','_mp_rbSN.csv')
    output_file = open(output_filename,'w')
    output_file.write(output_string)
    output_file.close()
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='input file')
    args = parser.parse_args()
    filename = args.i
    s = open(filename,'r')
    output_lines_list = []
    lines = s.readlines()
    data_lines = []
    for line in lines:
        data = line.split(',')
        if len(data) > 3:
            data_lines.append(line)
    print(len(data_lines))
    gff_lines = arab.open_gff()
    print('gff')
    pool_manager(data_lines,filename)
