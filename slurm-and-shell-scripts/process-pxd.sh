#!/bin/bash    

# load necesary modules
module load aspera
module load pgbtools

# initiates variables
D=0
M=0
F=0
pathtoparams=0
pathtofraggerslurms=0
pathtopythonqc=0
pathtodatabase=0
pxdlist=0

#
# This bit sets options and assigns arguments to 
# variables to be used later on
while getopts ':l:d:b:q:s:p:DMF' OPTION; do
  case "$OPTION" in
    l)
      pxdlist=$(realpath "$OPTARG")
      # check the pxd list exists
      if [ ! -f ${pxdlist} ]
      then
         echo "PXD list not found! Exiting..."
         exit 1
      fi
      if awk -F "\t" 'NF !=2' ${pxdlist}
      then
          printf "%-22s : %-40s %c\n" "PXD list" "${pxdlist}"
      else
          echo "PXD list is not a two column tab delimited list!"
          exit 1
      fi
      ;;
    
    b)
      pathtodatabase=$(realpath "$OPTARG")
      # check the protein database exists
      if [ ! -f ${pathtodatabase} ]
      then
          echo "Protein database not found! Exiting..."
          exit 1
      fi
  
      # check the protein database is a fasta
      if  head -1  ${pathtodatabase} | grep -q ">"
      then
          printf "%-22s : %-40s %c\n" "Protein database" "${pathtodatabase}" 
      else
          echo "Protein database not a fasta! Exiting..."
          exit 1
      fi
      ;;
    
    q)
      pathtopythonqc=$(realpath "$OPTARG")
      if [ -d ${pathtopythonqc} ] ; 
      then 
          printf "%-22s : %-40s %c\n" "QC scripts directory" "${pathtopythonqc}"
      else 
          echo "Missing QC directory Exiting..." 
          exit 1 
      fi
      
      ;;

    s)
      pathtofraggerslurms=$(realpath "$OPTARG")
      if [ -d ${pathtofraggerslurms} ] ; 
      then 
          printf "%-22s : %-40s %c\n" "MSfragger slurm script" "${pathtofraggerslurms}"
      else 
          echo "Missing MSfragger script directory Exiting..." 
          exit 1 
      fi
      ;;
    
    p)
      pathtoparams=$(realpath "$OPTARG")
      if [ -d ${pathtoparams} ] ; 
      then 
          printf "%-22s : %-40s %c\n" "MSfragger parameters" "${pathtofraggerslurms}"
      else 
          echo "Missing MSfragger parameters directory Exiting..." 
          exit 1 
      fi
      printf "%-22s : %-40s %c\n" "Parameters directory" "${pathtoparams}"
      ;;
    
    D)
      D=1
      ;;
      
    M)
      M=1
      ;;
      
    F)
      F=1
      ;;
    
    ?)
      echo "Unrecognised flag used."
      echo
      echo "Usage: $(basename $0) [-l pxd-list.txt] [-b database.fasta] [-q python-qc-directory] [-s slurm-script-directory] [-p fragger-params-directory] [-D] [-M] [-F]"
      echo
      echo "To run steps separately add one of the following flags:"
      echo "-D  run download PXD step"
      echo "-M  run convert raw to mzML step"
      echo "-F  run MSfragger and postprocessing step"
      exit 1
      ;;
  esac
done

echo

#echo $D
#echo $M
#echo $F
#echo $pathtoparams
#echo $pathtofraggerslurms
#echo $pathtopythonqc
#echo $pathtodatabase
#echo $pxdlist


# if the download option [-D] is given, execute this code:
if [ $D -eq 1 ]
then
tput setaf 3; echo "PXDs will be downloaded"; tput sgr0
for type in $(cut -f2 ${pxdlist} | sort -u)
      do 
          for i in $(grep -w ${type} ${pxdlist} | cut -f1)
          do
          tput setab 6; tput setaf 0; echo " ------ PROCESSING PXD:  ${i} ------ "; tput sgr 0
          echo

          # make a directory
          #mkdir ${i}
          #cd ${i} 
 
          # download the pride data
          tput setaf 6; echo "------ START of file download for ${i} ------"; tput sgr 0
	  echo
          pridepy download-all-public-raw-files -a ${i} -o ${i} -p aspera
          tput setaf 2; echo "------ END of file download for ${i} ------"; tput sgr 0
          echo

          n=$(ls -1q ${i}/*.raw 2> /dev/null | wc -l)
          if [ ${n} -eq 0 ]
          then
              tput setab 1; echo "No files for ${i} downloaded..." ; tput sgr 0
              echo
          fi
         #cd ../
         done
      done

fi


# if the modify option [-M] is given, execute this code:
if [ $M -eq 1 ]
then
tput setaf 3; echo "Raw files will be converted to mzML"; tput sgr0
    for type in $(cut -f2 ${pxdlist} | sort -u)
    do 
        for i in $(grep -w ${type} ${pxdlist} | cut -f1)
        do
        # navigate to PXD directory
        cd ${i}
        # get a list of raw files and the number of samples
        ls *.raw > samples.txt 2> /dev/null
        n=$(ls -1q *.raw 2> /dev/null | wc -l)
        if [ ${n} -eq 0 ]
        then
            tput setab 1; echo "No raw files for ${i} downloaded..." ; tput sgr 0 
            echo
        
        else
            # run mono in array to convert raw files to mzML
            tput setaf 6; echo "------ START of file conversion for ${n} raw files in ${i} ------"; tput sgr 0
            sbatch --array=1-${n} -c 2 -J ${i} --wait ${pathtofraggerslurms}/convert_mzML_array.sh 
            wait
            tput setaf 2; echo "------ END of file conversion for ${i} ------"; tput sgr0
            echo
        fi
        done
        # tidy up slurm outputs
        rm *.out
        rm *.err
        # navigate to parent directory
        cd ../
    done
fi


# if the MSfragger option [-F] is given, execute this code:

if [ $F -eq 1 ]
then
tput setaf 3; echo "MSfragger and quality checks will be run"; tput sgr0
      # set database in params files
      # database_name parameter should be on the first line of the params file
      for file in $(ls ${pathtoparams}/*)
      do
      p=$(echo "database_name = ${pathtodatabase}")
      sed -i -r "s|^(database_name =).*|${p}|" ${file}
      done
      
      # run fragger and qc scripts
      for type in $(cut -f2 ${pxdlist} | sort -u)
      do 
          for i in $(grep -w ${type} ${pxdlist} | cut -f1)
          do
          tput setab 6; tput setaf 0; echo " ------ PROCESSING PXD:  ${i} ------ "; tput sgr 0
          # navigate to directory
          cd ${i}
          n=$(ls -1q *.mzML 2> /dev/null | wc -l)
          if [ ${n} -eq 0 ]
          then
            tput setab 1; echo "No mzML files present for ${i}..." ; tput sgr 0 
            
          else
            tput setaf 6; echo "------ START of running MSfragger slurm for files in ${i} ------"; tput sgr 0
            echo "The type is set to ${type}"
            sbatch --array=1-1 ${pathtofraggerslurms}/TPP-fragger-wholedir-postprocess-GENERAL.sh ${pathtoparams}/fragger_${type}.params ${pathtopythonqc}
            tput setaf 3; echo "------ MSfragger slurm job for ${i} initiated ------"; tput sgr0
            echo
          fi

          cd ../
          done
      done

fi


## DEFAULT SCRIPT START if D, M, or F aren't given

if  [ $D -eq 0 -a $M -eq 0 -a  $F -eq 0  ]
then
    if [ "${pathtoparams}" == 0 ] || [ "${pathtofraggerslurms}" == 0 ] || [ "${pathtopythonqc}" == 0 ] || [ "${pathtodatabase}" == 0 ] || [ "${pxdlist}" == 0 ]
    then
        echo 'Missing one or more files and directories'
        echo
        echo "Usage: $(basename $0) [-l pxd-list.txt] [-b database.fasta] [-q python-qc-directory] [-s slurm-script-directory] [-p fragger-params-directory] [-D] [-M] [-F]"
        echo
        echo "To run steps separately add one of the following flags:"
        echo "-D  run download PXD step"
        echo "-M  run convert raw to mzML step"
        echo "-F  run MSfragger and postprocessing step"
        exit 1
    fi

    tput setaf 3; echo "Running default download, raw file conversion, and MSfragger processing"; tput sgr0 

    # set database in params files
    # database_name parameter should be on the first line of the params file
    for file in $(ls ${pathtoparams}/*)
    do
        p=$(echo "database_name = ${pathtodatabase}")
        sed -i -r "s|^(database_name =).*|${p}|" ${file}
    done

    # get the list of types
    for type in $(cut -f2 ${pxdlist} | sort -u)
    do 
        for i in $(grep -w ${type} ${pxdlist} | cut -f1)
        do
        tput setab 6; tput setaf 0; echo " ------ PROCESSING PXD:  ${i} ------ "; tput sgr 0
        
        # make a directory
        #mkdir ${i}
        #cd ${i}

        # download the pride data
        echo
        tput setaf 6; echo "------ START of file download for ${i} ------"; tput sgr 0; echo
        pridepy download-all-public-raw-files -a ${i} -o ${i} -p aspera
        tput setaf 2; echo "------ END of file download for ${i} ------"; tput sgr 0y

        echo

        # get a list of raw files and the number of samples
        cd ${i}
	ls *.raw > samples.txt 2> /dev/null
        n=$(ls -1q *.raw 2> /dev/null | wc -l )
        if [ ${n} -eq 0 ]
        then
            tput setab 1; echo "No raw files for ${i} downloaded..." ; tput sgr 0
            echo
        else
            # run mono in array to convert raw files to mzML
            tput setaf 6; echo "------ START of file conversion for ${n} raw files in ${i} ------"; tput sgr 0
            sbatch --array=1-${n} -c 2 -J ${i} --wait ${pathtofraggerslurms}/convert_mzML_array.sh 
            wait
            tput setaf 2; echo "------ END of file conversion for ${i} ------"; tput sgr0
            echo

            # run fragger and qc scripts
            tput setaf 6; echo "------ START of running MSfragger slurm for files in ${i} ------"; tput sgr 0
            echo "The type is set to ${type}"
            sbatch --array=1-1 ${pathtofraggerslurms}/TPP-fragger-wholedir-postprocess-GENERAL.sh ${pathtoparams}/fragger_${type}.params ${pathtopythonqc}
            tput setaf 3; echo "------ MSfragger slurm job for ${i} initiated ------"; tput sgr0
            echo
            cd ../
        fi
        done
    done
fi


