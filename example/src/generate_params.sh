#!/opt/homebrew/bin/bash
# Sergi Ortiz 6 Apr 2025
# Adapted from MolBioMed GitHub
# preprocess the extracted pdb such that everything works

for f in *.pdb; do

    # define variables
    echo "$f"
    echo "${f%.pdb}"
    prmtop_file="${f%.pdb}_solv.prmtop"
    inpcrd_file="${f%.pdb}_solv.inpcrd"

    # echo variables
    echo "${prmtop_file}"

    mkdir ${f%.pdb}
    mv $f ${f%.pdb}
    cd ${f%.pdb}
    cat $f | grep -v EPW > ${f%.pdb}_TIP3P.pdb.0
    cat ${f%.pdb}_TIP3P.pdb.0 | grep -v CONECT > ${f%.pdb}_TIP3P.pdb.1
    sed "s/359 X/HD1 X/g" ${f%.pdb}_TIP3P.pdb.1 > ${f%.pdb}_TIP3P.pdb.2
    sed "s/364 X/HD2 X/g" ${f%.pdb}_TIP3P.pdb.2 > ${f%.pdb}_TIP3P.pdb.3
    sed "s/539 X/HD3 X/g" ${f%.pdb}_TIP3P.pdb.3 > ${f%.pdb}_TIP3P.pdb.4
    sed "s/543 X/HE1 X/g" ${f%.pdb}_TIP3P.pdb.4 > ${f%.pdb}_TIP3P.pdb.5
    sed "s/661 X/IE1 X/g" ${f%.pdb}_TIP3P.pdb.5 > ${f%.pdb}_TIP3P.pdb
    #pdb4amber -i ${f%.pdb}_TIP3P_final.pdb -o ${f%.pdb}_TIP3P.pdb

    rm ${f%.pdb}_TIP3P.pdb.*
    cp -r ../tleap .
    cd tleap
    mv ../${f%.pdb}_TIP3P.pdb .

    sed "s/{frame}/${f%.pdb}/g" tleap.template > tleap.inp

    tleap -f tleap.inp | tee tleap.out

    echo "using `which python`"

    # flags
    # -cs -csc 665 -csr 17              crop solvent around residue 665 with a radius of 17 A. 
    # -al -alc 10612 -alr 15            create active atom list around atom number 10612 with a radius of 15 A.
    # -ra -ram N N N N O O O FE -rac    crop and modify prmtop topology so that ChemShell understands. Checks if any atom is unknown by chemshell
    #               path to python script               prmtop topology file            inpcrd coordinate file
    python ../../amber_to_chemshell_modeller.py -p ${f%.pdb}_solv.prmtop -c ${f%.pdb}_solv.inpcrd -cs -csc 665 -csr 17 -al -alc <CATOM> -alr 15 -ra -ram N N N N O O O FE -rac

    cd ../../
done
