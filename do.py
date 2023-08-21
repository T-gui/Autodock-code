import sys
import os
import pymol
from pymol import cmd
import itertools
NATURAL_AMINO_ACIDS = ["ARG", "ALA", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LED", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
def residue_mutation(protein,chain,mutation_positions, output_dir):
    #t该函数目的在于产生result文件夹以及产生组合数列表
    total_positions = len(mutation_positions)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for r in range(1, total_positions+1):
        combinations = itertools.combinations(mutation_positions, r)
        #产生组合数列表，并对每一组组合数单独计算
        for combination in combinations:
            newdir = os.path.join(output_dir,str(combination))
            divide_mutation(protein,chain,combination,newdir)
    


def divide_mutation(protein,chain,combination,output_dir):
    #该方法为对每一个组合数序列进行依次突变，每一次突变都会产生一个文件夹，下一次突变将从上一个文件夹中选取蛋白质，以保证突变的完整性
    position = 0
    filename = output_dir
    positions = len(combination)
    while position < positions :
        #按照组合数中的数字依次进行突变
        position_name = '%s'% position
        last_position = '%s' % (position-1)
        filenames= os.path.join(filename,position_name)
        if not os.path.exists(filenames):
            os.makedirs(filenames)
        last_file = os.path.join(filename,last_position)
        if not os.path.exists(last_file):
            #如果该次突变是第一次突变，则直接对该位点进行20次突变
            for aa in NATURAL_AMINO_ACIDS:
                cmd.load(protein,'pro_mutation')
                cmd.wizard("mutagenesis")
                cmd.do("refresh_wizard")
                cmd.get_wizard().set_mode(aa)
                cmd.get_wizard().do_select("/pro_mutation//%s/%d" % (chain, combination[position]))
                cmd.frame(1)
                cmd.get_wizard().apply()
                output_filename = os.path.join(filenames, "%s_%d_%s_%s" % (chain,combination[position], aa, protein))
                cmd.save(output_filename, 'pro_mutation')
                cmd.delete("all")
        if os.path.exists(last_file):
            #如果该次突变不是第一次突变，则从上一次突变文件夹中选取突变目标，进行突变
            file_names=os.listdir(last_file)
            for aa in file_names:
                old_position = os.path.join(last_file,aa)
                mutation(old_position,filenames,combination[position],chain,protein)
            #突变完毕，删除上一个文件夹以节省空间
            for root, dirs, files in os.walk(last_file, topdown=False):
                for name in files:
                    os.remove(os.path.join(root, name))
                for name in dirs:
                    os.rmdir(os.path.join(root, name))
            os.removedirs(last_file)
        position = position+1
def mutation(oldfile,newfile,index,chain,protein):
    #该方法抽象了突变过程，引用该方法即可完成一次突变
    for aa in NATURAL_AMINO_ACIDS :
        cmd.load(oldfile,'pro_mutation')
        cmd.wizard("mutagenesis")
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode(aa)
        cmd.get_wizard().do_select("/pro_mutation//%s/%d" % (chain, index))
        cmd.frame(1)
        cmd.get_wizard().apply()
        proteins = os.path.basename(oldfile)
        newname = proteins.replace(protein,"")
        combination_str = newname+"%d_%s_%s" % (index, aa, protein)
        output_filename = os.path.join(newfile,combination_str)
        cmd.save(output_filename, 'pro_mutation')
        cmd.delete("all")
        
def main():
    protein = str(sys.argv[1])
    chain = str(sys.argv[2])
    mutation_positions = [int(p) for p in sys.argv[3].split(',')] 
    output_dir = "result"  # Name of the output directory
    residue_mutation(protein,chain,mutation_positions,output_dir)

if __name__=="__main__":
    main()