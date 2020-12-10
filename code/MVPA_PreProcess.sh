#!/bin/sh
#
#$ -cwd
#$ -o output/o_PreP_CoRe_001.txt
#$ -e output/e_PreP_CoRe_001.txt
#$ -m bae

# export PATH=$PATH:/cerebro/cerebro1/dataset/core_mvpa/code/toolbox/afni
#
# export FSLDIR=/cerebro/cerebro1/dataset/core_mvpa/code/toolbox/fsl
# source /cerebro/cerebro1/dataset/core_mvpa/code/toolbox/fsl/etc/fslconf/fsl.sh

export PATH=$PATH:/export01/local/afni 
#
export FSLDIR=/export01/local/fsl
source ${FSLDIR}/etc/fslconf/fsl.sh

usage() {
cat<<EOF

==================== { MVPA_PreProcess.sh } ========================================
||
||  [Purpose]: Loop through each sleep (in Day 1 & 2 & 3) and all other task/resting sessions of each subjects to do the preprocessing for CoRe dataset:
||
||             -----[   1.  Motion corretion or realignment;  ]
||             -----[   2.  Despiking;  ]
||             -----[   3.  Coregistration to sleep(@D1), segementation and inverse normalization to prepare native masks for all ROIs and CSF, WM;  ]
||                          ***** <  optional: 3dwarp deoblique the reference EPI image, thereafter all other sessions datasets through coregistration >
||             -----[   4.  Noise removal: Trends(1/2) + Motion(6*4) + CSF/WM(2*4) = 33/34 parameters;  ]
||                          ***** <  caution: only sleep session has quadratic trend to be removed  >
||             -----[   5.  Resampling with/without Normalization;  ]
||                          ***** <  optional: resample the data  >
||             -----[   6.  Smoothing;  ]
||             -----[   7.  Z-scoring;  ]
||
||  [Usage]: bash MVPA_PreProcess.sh

||           < specify the vairable array of "sesslist", "sessshort", "daylist", "subj" below >
||           < modify the parameters below >
||
||  [Notice]: a)
||
===========================================================================================
Author: Shuo, Chen (shuo.chen2@mail.mcgill.ca)
  Date: 2020.0102 (bash shell)
      @ JD lab, BIC of MNI, McGill University, Canada

EOF
exit
}

if [ "$1" == '-help' ]; then
    usage
    exit 1
fi


# $ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ==================  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$      MAIN CODES      $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# $ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  ==================  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#=================== [ specify subjects and the parent directory ] ====================#

datapath="/cerebro/cerebro1/dataset/core_mvpa/"

declare -a subjlist=("CoRe_001") # CoRe_001 CoRe_023 CoRe_063 CoRe_079 CoRe_107 CoRe_128 CoRe_082 CoRe_100 CoRe_087 CoRe_192 CoRe_195 CoRe_220 CoRe_223 CoRe_235 CoRe_237 CoRe_267 CoRe_268 CoRe_296 CoRe_011 CoRe_022 CoRe_050 CoRe_054 CoRe_094 CoRe_102 CoRe_155 CoRe_162


#=================== [ specify file and data structure when doing the fine tuning] ====================#

AutoSearch=('n') # y/n

# ### all sessions for a given subject.
# # declare -a sesslist=("12-sleep" "11-sleep" "03-rest" "04-TSeq" "05-rest" "06-rest" "07-Retest_TSeq" "04-rest" "05-Retest_TSeq" "06-IntSeq" "07-rest" "11-sleep" "03-rest" "04-TSeq" "05-rest" "06-IntSeq" "07-MVPA1" "08-MVPA2")
# # declare -a sessshort=("sleep" "sleep" "rest" "TSeq" "rest" "rest" "Retest_TSeq" "rest" "Retest_TSeq" "IntSeq" "rest" "sleep" "rest" "TSeq" "rest" "IntSeq" "MVPA1" "MVPA2")
# # declare -a daylist=("D1" "D1" "D1" "D1" "D1" "D1" "D1" "D2" "D2" "D2" "D2" "D2" "D3" "D3" "D3" "D3" "D3" "D3")
sleep_ref="11-sleep" # "15-sleep"

# # declare -a sesslist=("04-rest" "05-Retest_TSeq" "06-IntSeq" "07-rest" "11-sleep")
# # declare -a sessshort=("rest" "Retest_TSeq" "IntSeq" "rest" "sleep")
# # declare -a daylist=("D2" "D2" "D2" "D2" "D2")
#
declare -a sesslist=("11-sleep" ) #
declare -a sessshort=("sleep" ) #
declare -a daylist=("D1")

#=================== [ modify handles and parameters] ====================#
# A same preprocessing pipeline be applied to all resting sessions

Loops2go=('2') # all, 1, 2, 3, chk : all = go all loops; otherwise each number specify single loop to go; chk = checking the result for coregistration and normalization

Motion_Outlier=('0.3') # threshold for outliers, calculated by afni enorm
Deoblique=('n'); # deoblique the reference EPI image from sleep session, but later sleep session need to be coregistrated on it again
Pause_Dspike=('n') # y = yes(do not conduct despiking AGAIN); skip despiking in loop 2; mostly in order to fine tune the parameters for the preprocessing of anatomical image
Shrink_Box=('y'); # y/n : shrink the box size for every EPI sessions after despiking, to the same size as the cropped coreg_base image's
Pause_Coreg=('n') # y = yes(do not conduct coregistration AGAIN); skip coregistration in loop 2; mostly in order to fine tune the parameters for the preprocessing of anatomical image

T1_resolu=('1mm'); # 1mm or 2mm. Options for anatomical image resolution and the corresponding MNI template resolution
BET_thre=('0.35') # threshold for scull removal; 0.32 preferred
Seg_prob=('.95') # CAUTION : use the format .xx not 0.xx; for segmentation
NormHandle=('n') # l/n : linear or nonlinear normalization
Masklist=('csf' 'wm'); # Masks used for physiological noise
WM_option=('invnorm'); # invnorm/seg : invnorm = using inverse-normalized mask; seg = using natively segmented mask
CSF_option=('seg'); # as above row
CSF_noborder=('y'); # given CSF_option=seg, considering to discard the border voxels
InvNmask_thre=('.95') # CAUTION : use the format .xx not 0.xx; for inverse normalization

CMthreshold=('6') # mm ; if center distance > threshold, then do the center matching before the coregistration between EPI sessions
Smth_kernel=('6') # mm
Resample_Norm=('rn'); # r = only resampling at native space; rn = resampling+normalization; ne = neither
Final_voxelsize=('3') # o = keep original as EPI dataset; n = specify a size number mm^3; Resampling

ReCheck=('n'); # y = delete the checking folders, then recreate; n = not
# Match_origin=('y') # y/n : yes or no to match center before coregistration between epi and anat images



#=========================== [ lets LOOP it up ] ============================#



for subj in "${subjlist[@]}"; # basically, 3 loops cascade sequentially
do

    if [ "${AutoSearch}" == "y" ] ;
    then

      # to obtain the file structure for each subject
      sesslist=()
      sessshort=()
      daylist=()
      sleep_ref=''
      sleep_vol_num=0
      declare -a all_3_days=("D1" "D2" "D3")

      for day in "${all_3_days[@]}" ;
      do
        cd ${datapath}/MRI_nii/${subj}/$day

        for ss in */ ;
        do
          ss_cut=$(echo $ss | cut -d '-' -f 2)

          if [ "${ss_cut}" == "rest/" ] || [ "${ss_cut}" == "sleep/" ] || [ "${ss_cut}" == "TSeq/" ] || [ "${ss_cut}" == "Retest_TSeq/" ] || [ "${ss_cut}" == "IntSeq/" ] || [ `echo $ss | rev | cut -c3-6 | rev` == "MVPA" ];
          then
            # echo "${ss}"

            sesslist+=($(echo $ss | cut -d '/' -f 1))
            sessshort+=($(echo $(echo $ss | cut -d '/' -f 1) | cut -d '-' -f 2))
            daylist+=($day)

          fi

          if [ "${day}" == "D1" ] && [ "${ss_cut}" == "sleep/" ] ;
          then
            sleep_vol_temp=`3dinfo -nv ${ss}*sleep.nii`
            # echo "${sleep_vol_temp}"

            if [ `echo "${sleep_vol_temp} > ${sleep_vol_num}" | bc` -eq 1 ]; # choose the longest sleep scan to have the EPI_ref image
            then
              sleep_vol_num=${sleep_vol_temp}
              sleep_ref=($(echo $ss | cut -d '/' -f 1))
            fi
          fi

        done
      done

      # echo "${sesslist[@]}"
      # echo "${sessshort[@]}"
      # echo "${daylist[@]}"
      # echo "${sleep_ref}"

    fi


      # get length of an sesslist
      sesslength=${#sesslist[@]}
      # sesslength=1

      workdir="${datapath}/MRI_PreP"

      if [ ! -d "${datapath}/MRI_PreP/${subj}" ]; then
        echo
        echo "|| ----- First time handle, copying files into ${datapath}/MRI_PreP/${subj} -----"
        echo
        cp -R ${datapath}/MRI_nii/${subj} ${datapath}/MRI_PreP
      fi



  #   # ------------------------------------------------------------------------------------------------------------
  #   # ^w^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ <<  loop 1  >> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^o^

      if [ "${Loops2go}" == "all" ] || [ "${Loops2go}" == "1" ]; then
        for (( ss=0; ss<${sesslength}; ss++ ));      # because start from ss=0, and in bash index 0=1st item in the array
        # for (( ss=1; ss<${sesslength}+1; ss++ ));  # then all should be replaced like this ${sesslist[$ss-1]}
        do

          cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}
          data="${subj}_${daylist[$ss]}_${sessshort[$ss]}"

          echo ' '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo '|| '
          echo "||    -----[   1. motion corretion or realignment;  ]-----"
          echo '|| '
          echo "||    current session under preprocessing is : [ '${sesslist[$ss]}' of ${subj} ]"
          echo "||                    the input file name is : [ ${data}.nii ]"
          echo '|| '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo ' '

          # find the volume/time-point with minimal voxel counts of outliers, to serve as reference base
          3dToutcount -automask -range -polort 1 "${data}.nii" > "${data}_outlier.txt" # other parameters: -fraction -legendre
          base=`cat ${data}_outlier.txt | perl -0777an -F"\n" -e '$i=0; $small=999999; map {/\s*(\d+)/; if ($small > $1) {$small = $1; $ind=$i}; $i++;} @F; print $ind'` # search for the minimal movement frame?

          # for later coregistration use
          if [ "${daylist[$ss]}" == "D1" ] && [ "${sesslist[$ss]}" == "${sleep_ref}" ];
          then

            echo
            echo "|| ----- the base reference is volume ${base} -----"
            echo

            rm ${subj}_coreg_base.nii.gz ${subj}_coreg_base_*.nii.gz
            3dcalc -a "${data}.nii[$base]" -expr "a" -prefix "${subj}_coreg_base_beforeAutoboxBackup.nii.gz" # because averaging the long sleep session to have mean ref_image, instead, may dilate the border of brain
            # 3dcalc -a "${data}.nii[$base]" -expr "a" -prefix "${subj}_coreg_base.nii.gz" # because averaging the long sleep session to have mean ref_image, instead, may dilate the border of brain

            EPIbase="${subj}_coreg_base"

            # editing the box for EPI ref, so that after coregistration the T1 image will inherit this box to speed up the nonlinear normalization
            3dAutobox -input "${EPIbase}_beforeAutoboxBackup.nii.gz" -prefix "${EPIbase}_temp1.nii.gz"  # squeeze the box; without the -noclust, the box would be really tight
            3dZeropad -I 2 -S 2 -A 3 -P 3 -L 3 -R 3 -prefix "${EPIbase}_temp2.nii.gz" "${EPIbase}_temp1.nii.gz"
            3dZeropad -prefix "${EPIbase}.nii.gz" -master "${EPIbase}_temp2.nii.gz" "${EPIbase}_beforeAutoboxBackup.nii.gz"
            rm ${EPIbase}_temp*.nii.gz


            if [ "${Deoblique}" == "y" ]; then # all other sessions and anatomical will be coregistered upright with the ref_plumb image later
              3dZeropad -RL 100 -AP 100 -IS 80 -prefix "${EPIbase}_pad.nii.gz" "${EPIbase}.nii.gz" # dilate the box, so that after 3dWarp the image is still in the box
              3dWarp -deoblique -prefix ${EPIbase}_plumb_temp1.nii.gz -gridset ${EPIbase}_pad.nii.gz ${EPIbase}.nii.gz # make sure no resampling(but you can do the resampling here by -newgrid), and keep brain intact inside the box
              3dAutobox -input ${EPIbase}_plumb_temp1.nii.gz -prefix ${EPIbase}_plumb_temp2.nii.gz -noclust # squeeze the box
              3dZeropad -RL 64 -AP 64 -IS 40 -prefix ${EPIbase}_plumb.nii.gz ${EPIbase}_plumb_temp2.nii.gz # recover the box
              # rm ${EPIbase}_plumb_temp*.nii.gz ${EPIbase}_pad.nii.gz
              EPIbase="${EPIbase}_plumb"
            fi
          fi

          # motion correction : coregistrate all volumes to the base volume; ${data}m, M = motion
          3dvolreg -base ${base} -twopass -1Dfile "${data}_motion.txt" -prefix "${data}_M.nii.gz" "${data}.nii" # other parameters: -verbose -zpad 1 -heptic -1Dmatrix_save mat.r01.1D

          # create a temporal mask to mark down huge instant (delta/relative) movement volumes; motion censoring/scrubbing
          1d_tool.py -infile "${data}_motion.txt" -set_nruns 1 -show_censor_count -censor_prev_TR -censor_motion ${Motion_Outlier} ${data} # -censor_next_TR -censor_prev_TR
            # try stricter threshold as 0.2 for both sleep and task sessions. More scrubbing no harm, as long we have enough samples for both training and prediction
            # https://afni.nimh.nih.gov/afni/community/board/read.php?1,161451,161463#msg-161463
            # https://afni.nimh.nih.gov/afni/community/board/read.php?1,143809,144115
            # <if scrubbed out beyond 50% volumes, discard the session>

            # in case other indice for calculating to-be-scrubbed volumes are required
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_DV_ConfM_${subj}_${session}.txt --dvars -s Mo_DV_raw_${subj}_${session}.txt -p Mo_DV_Plot_${subj}_${session}  #  based on intensity differences within the realigned timeseries; In the case that the motion correction is not accurate, then using motion correction parameters (rotation angles and translations in mm) is a poor way to estimate the outliers.
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_FD_ConfM_${subj}_${session}.txt --fd -s Mo_FD_raw_${subj}_${session}.txt -p Mo_FD_Plot_${subj}_${session} # Power
              # fsl_motion_outliers -i "${data}.nii.gz" -o Mo_FDrms_ConfM_${subj}_${session}.txt --fdrms -s Mo_FDrms_raw_${subj}_${session}.txt -p Mo_FDrms_Plot_${subj}_${session} # Jenkinson
              # <compare fsl --refrms --fdrms and afni enorm .3>

          # rm "${data}_outlier.txt"
          unset base
        done

      # else          ## not implementing loop 1 AGAIN, for the convenience to pass down settings

      fi


      #         ---- << code between loop >> ----

      ## [ SET UP ] : some vairables here between loop1 and loop2, in case there is no sleep session included in a RUN

      # base_path="${workdir}/${subj}/D1/`ls -d *sleep`" # damn, the sleep session number sometimes changes; D1 for sure; e.g. CoRe_082 become 12-sleep
      base_path="${workdir}/${subj}/D1/${sleep_ref}" # damn, some sbj's sleep has two sessions, need to be specified; e.g. CoRe_001, has both 11-sleep and 12-sleep

      EPIbase="${subj}_coreg_base"
      T1_hr="${subj}_anat_${T1_resolu}"
      if [ "${Deoblique}" == "y" ]; then # all other sessions and anatomical will be coregistered upright with the ref_plumb image later
        EPIbase="${EPIbase}_plumb"
      fi
      cd ${workdir}/${subj}/D1

      #         ---- << code between loop >> ----


      # ------------------------------------------------------------------------------------------------------------
      # ^w^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ <<  loop 2  >> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^o^
      if [ "${Loops2go}" == "all" ] || [ "${Loops2go}" == "2" ]; then

        for (( ss=0; ss<${sesslength}; ss++ ));      # because start from ss=0, and in bash index 0=1st item in the array
        do

          cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}
          data="${subj}_${daylist[$ss]}_${sessshort[$ss]}_M"

          echo ' '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo '|| '
          echo "||    -----[   2. Despiking ;  ]-----"
          echo "||    -----[   3. Coregistration to sleep(@D1), preparing CSF,  WM and ROIs masks ;  ]-----"
          echo '|| '
          echo "||    current session under preprocessing is : [ '${sesslist[$ss]}' of ${subj} ]"
          echo "||                    the input file name is : [ ${data}.nii.gz ]"
          echo '|| '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo ' '

          ## Despiking first
          if [ "${Pause_Dspike}" == "y" ]; then
            echo
            echo "|| ----- Pausing the despiking step, no need for repetition -----"
            echo
          else
            rm "${data}Ds.nii.gz" "${data}d_spike.nii.gz"
            3dDespike -ssave "${data}d_spike.nii.gz" -prefix "${data}Ds.nii.gz" "${data}.nii.gz"
          fi

          data="${data}Ds"


          # shrink the box for all EPI sessions as the same as the coreg base image.
          if [ "${Shrink_Box}" == "y" ]; then
            rm "${data}_bSBbackup.nii.gz"
            3dcopy "${data}.nii.gz" "${data}_bSBbackup.nii.gz"

            if [ "${daylist[$ss]}" == "D1" ] && [ "${sesslist[$ss]}" == "${sleep_ref}" ];
            then
              3dZeropad -prefix "${data}_smallbox.nii.gz" -master "${EPIbase}.nii.gz" "${data}.nii.gz"
            else
              cp "${base_path}/${EPIbase}.nii.gz" "${EPIbase}.nii.gz"
              3dZeropad -prefix "${data}_smallbox.nii.gz" -RL `3dinfo -ni *base.nii.gz` -AP `3dinfo -nj *base.nii.gz` -IS `3dinfo -nk *base.nii.gz` "${data}.nii.gz" # -master still could fail to match the slices number
            fi

            rm "${data}.nii.gz"
            3dcopy "${data}_smallbox.nii.gz" "${data}.nii.gz"
            rm "${data}_smallbox.nii.gz"
          fi


          # because we only need to work on anatomical image ONCE, take the convenience here,
          # whenever the session equal to the EPI reference session, finish the preprocess anatomical image as well
          if [ "${daylist[$ss]}" == "D1" ] && [ "${sesslist[$ss]}" == "${sleep_ref}" ];
          then


            # < steps inside this module >
            # 1. Anat scull removal; 2. coreg to sleep base; 3. seg in native space to have CSF mask; 4. norm to mni, then apply inverse xform for having WM mask

            ## 2do: fsl fnirt seems working; compare to afni 3dQwarp (newest wu code) for nonlinear normalization

            cd ${workdir}/${subj}/anat_HR
            rm ${subj}_anat_*mm.nii.gz

            # 3dcopy ${workdir}/${subj}/anat/${subj}_anat.nii ${subj}_anat_2mm.nii.gz # commented because it is hard to extract brain from 2mm.T1 using BET but surprisingly working with 3dSkullStrip
            3dcopy ${subj}_anat_HR.nii ${subj}_anat_1mm.nii.gz


            # 1. Anat scull removal;

            rm *+Bet*.nii.gz

            bet "${subj}_anat_1mm.nii.gz" "${subj}_anat_1mm+Bet.nii.gz" -f ${BET_thre} -g 0 -B -R

            # bet "${T1_hr}.nii.gz" "${T1_hr}+Bet.nii.gz" -R -f ${BET_thre} -g 0 -B #-F # used to try scull removal directly on 2mm anatomical images
            # bet "${EPIbase}.nii.gz" "${EPIbase}_bet.nii.gz" -R -f 0.2 -g 0 # bet the epi base image
            # 3dSkullStrip -input CoRe_011_anat_2mm.nii.gz -push_to_edge -touchup -ld 30 -prefix "test+Bet_afni.nii.gz"
            # [special for sbj CoRe_128]
              # 3dSkullStrip -input "${subj}_anat_1mm.nii.gz" -shrink_fac 0.45 -shrink_fac_bot_lim 0.5 -touchup -touchup -push_to_edge -max_inter_iter 10 -ld 20 -prefix "${subj}_anat_1mm+Bet.nii.gz" # -touchup -push_to_edge -touchup -shrink_fac 0.85 -shrink_fac_bot_lim 0.7 -no_touchup -no_push_to_edge -max_inter_iter 10 -interactive -fac 1.5

            if [ "${T1_resolu}" == "2mm" ]; then
              echo
              echo "|| ----- downsampling for a 2 mm^3 scull-removed anatomical image for ${subj} -----"
              echo
              3dresample -dxyz 2 2 2 -prefix "${subj}_anat_2mm+Bet.nii.gz" -input "${subj}_anat_1mm+Bet.nii.gz"
              # [special for sbj CoRe_128]
              #   3dcopy "${subj}_anat_1mm+Bet.nii.gz" "${subj}_anat_2mm+Bet.nii.gz" # squeeze the box for later nonlinear normalization

            fi

            T1_hr="${T1_hr}+Bet"

            # #   2do: go next step, see whether affect normalization

            3dcopy ${T1_hr}.nii.gz ${T1_hr}_beforeAutoboxBackup.nii.gz
            rm ${T1_hr}.nii.gz
            3dAutobox -input "${T1_hr}_beforeAutoboxBackup.nii.gz" -prefix "${T1_hr}.nii.gz" -noclust # squeeze the box for later nonlinear normalization

            # 2. Coreg to sleep base;
            rm *${EPIbase}*
            cp -f "${base_path}/${EPIbase}.nii.gz" "${EPIbase}.nii.gz"

            # 3dcopy ${EPIbase}.nii.gz ${EPIbase}_beforeAutoboxBackup.nii.gz
            # rm ${EPIbase}.nii.gz
            # 3dAutobox -input "${EPIbase}_beforeAutoboxBackup.nii.gz" -prefix "${EPIbase}.nii.gz" -noclust # squeeze the box

            T1_path=`pwd`


            CMDist="$(@Center_Distance -dset ${T1_hr}.nii.gz ${EPIbase}.nii.gz)"

            echo
            echo "|| ----- CM distance = ${CMDist} mm -----"
            echo

            rm *_epicoreg+*


            align_epi_anat.py -anat2epi -epi2anat -ginormous_move -anat "${T1_hr}.nii.gz" -epi "${EPIbase}.nii.gz" \
            -epi_base 0 -suffix "_epicoreg" -anat_has_skull no -epi_strip 3dAutomask -cost "lpc+ZZ" \
            -ex_mode quiet -volreg off -tshift off -save_resample -align_centers on -Allineate_opts '-weight_frac 1.0 -maxrot 6 -maxshf 10 -VERB -warp shr'


            rm *${T1_hr}Co*
            3dcalc -a "${T1_hr}_epicoreg+*.HEAD" -expr 'a' -prefix "${T1_hr}Co.nii.gz"
            # 3dresample -master "${base_path}/${data}.nii.gz" -prefix "${T1_hr}Co_rs.nii.gz" -inset "${T1_hr}Co.nii.gz"


            # 3dcalc -a "${T1_hr}_epicoreg+*.HEAD" -expr 'a' -prefix "${T1_hr}Co_temp.nii.gz"
            # 3dAutobox -input "${T1_hr}Co_temp.nii.gz" -prefix "${T1_hr}Co.nii.gz" -noclust # squeeze the box to save time during nonlinear normalization
            # rm "${T1_hr}Co_temp.nii.gz"

            rm *_epicoreg+*

            T1_hr="${T1_hr}Co"


            # 3. seg in native space to have CSF mask;

            fast --channels=1 --type=1 --class=3 -B ${T1_hr}.nii.gz #--out CoRe_011_bet

            # by default, pve_0=CSF
            rm ${subj}_mask_*_seg.nii.gz
            3dcalc -a ${T1_hr}_pve_0.nii.gz -expr "ispositive(a-${Seg_prob})" -prefix "${T1_hr}_pve_0_thr${Seg_prob}.nii.gz"
            3dmask_tool -input ${T1_hr}_pve_0_thr${Seg_prob}.nii.gz -prefix ${T1_hr}_pve_0_thr${Seg_prob}.erod1.nii.gz -dilate_input -1
            3dcalc -a ${T1_hr}_pve_0_thr${Seg_prob}.erod1.nii.gz -expr "a" -prefix ${subj}_mask_csf_seg.nii.gz

            # by default, pve_2=WM
            3dcalc -a ${T1_hr}_pve_2.nii.gz -expr "ispositive(a-${Seg_prob})" -prefix "${T1_hr}_pve_2_thr${Seg_prob}.nii.gz"
            3dmask_tool -input ${T1_hr}_pve_2_thr${Seg_prob}.nii.gz -prefix ${T1_hr}_pve_2_thr${Seg_prob}.erod1.nii.gz -dilate_input -1
            3dcalc -a ${T1_hr}_pve_2_thr${Seg_prob}.erod1.nii.gz -expr "a" -prefix ${subj}_mask_wm_seg.nii.gz


            # 4. norm to mni, then apply inverse xform for having WM mask (because some sub cortical structures like thalamus are always be mostly regarded as WM in seg masks )

            for physio_mask in "${Masklist[@]}";
            do

              rm -f mask_*${physio_mask}*.nii.gz *${physio_mask}*_invnorm.nii.gz
              cp -f ${datapath}/Analysis_cabinet/mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_${physio_mask}_tal_nlin_asym_09c_thr.95*.nii.gz mask_${physio_mask}%mni.nii.gz
              # thr.95*.gz : CSF use the un-eroded mni mask; WM use the eroded one. original files were added with .hidden suffix in the original folder

              image_to_move="mask_${physio_mask}%mni"



              # just need to be done once for the t1hr2mni normalization; thereafter all the masks will use the same inv-matrix
              if [  "${physio_mask}" == "${Masklist[0]}" ]; then

                echo
                echo '|| ----- Generating new transformation matrix for Inverse Transform... [FSL] -----'
                echo

                cp -f "${datapath}/Analysis_cabinet/mnit1_09c_${T1_resolu}.nii.gz" "mnit1_09c_${T1_resolu}.nii.gz"

                flirt -in "${T1_hr}.nii.gz" -ref "${datapath}/Analysis_cabinet/mnit1_09c_${T1_resolu}.nii.gz" -omat "t1hr2mni.mat" -out "${T1_hr}%mni.nii.gz"\
                    -bins 256 -cost corratio -searchrx 0 0 -searchry 0 0 -searchrz 0 0 -dof 12
                convert_xfm -omat "t1hr2mni_inv.mat" -inverse "t1hr2mni.mat"

                if [  "${NormHandle}" == "n" ]; then
                  rm t1hr2mni_NLwarp*
                  fnirt --in="${T1_hr}.nii.gz" --aff=t1hr2mni.mat --cout=t1hr2mni_NLwarp --config=T1_2_MNI152_2mm --iout="${T1_hr}%mniNL.nii.gz" \
                      --ref="${datapath}/Analysis_cabinet/mnit1_09c_${T1_resolu}.nii.gz" \
                      --refmask="${datapath}/Analysis_cabinet/mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_t1_tal_nlin_asym_09c_mask_${T1_resolu}.nii"
                  invwarp --ref=${T1_hr}.nii.gz --warp=t1hr2mni_NLwarp --out=t1hr2mni_NLwarp_inv
                fi
              fi


              # apply inv norm to physio_masks; either in a linear or nonlinear way
              if [  "${NormHandle}" == "l" ]; then
                flirt -ref "${T1_hr}.nii.gz" -in "${image_to_move}.nii.gz" -applyxfm -init "t1hr2mni_inv.mat"\
                    -schedule ${FSLDIR}/etc/flirtsch/xyztrans.sch -usesqform \
                    -omat "xform_temp.mat" -out "xform_temp.nii.gz"
                flirt -in "${image_to_move}.nii.gz" -ref "xform_temp.nii.gz" -noresampblur -noresample -noclamp \
                    -applyxfm -init  "xform_temp.mat" -out "mask_${physio_mask}%mni#native.nii.gz"

              elif [  "${NormHandle}" == "n" ]; then
                applywarp --ref="${T1_hr}.nii.gz" --in="${image_to_move}.nii.gz" \
                    --warp=t1hr2mni_NLwarp_inv --out="mask_${physio_mask}%mni#native.nii.gz" --interp=nn # "nearest neighbour" interpolation : non-integer values along the edges of the ROI's; suitable for single-vaule mask (1,0)
              fi


              # thresholding and renaming the resultant mask
              3dcalc -a "mask_${physio_mask}%mni#native.nii.gz" -expr "ispositive(a-0${InvNmask_thre})" -prefix "mask_${physio_mask}%mni#native.thr${InvNmask_thre}.nii.gz"
              3dcalc -a "mask_${physio_mask}%mni#native.thr${InvNmask_thre}.nii.gz" -expr "a" -prefix ${subj}_mask_${physio_mask}_invnorm.nii.gz

              # prepare a dilated mask for csf from segmentation, to eliminate csf voxels in the border of brain; no worry, you will have these transfrom files during the 1st run of this for-loop
              if [  "${CSF_option}" == "seg" ] && [  "${CSF_noborder}" == "y" ]; then
                cp -f ${datapath}/Analysis_cabinet/mni_icbm152_nlin_asym_09c_nifti/mni_icbm152_csf_tal_nlin_asym_09c_thr.95.dila1.nii mask4csf%mni.nii.gz
                image_to_move="mask4csf%mni"
                if [  "${NormHandle}" == "l" ]; then # apply inv norm to physio_masks; either in a linear or nonlinear way
                  flirt -ref "${T1_hr}.nii.gz" -in "${image_to_move}.nii.gz" -applyxfm -init "t1hr2mni_inv.mat"\
                      -schedule ${FSLDIR}/etc/flirtsch/xyztrans.sch -usesqform \
                      -omat "xform_temp.mat" -out "xform_temp.nii.gz"
                  flirt -in "${image_to_move}.nii.gz" -ref "xform_temp.nii.gz" -noresampblur -noresample -noclamp \
                      -applyxfm -init  "xform_temp.mat" -out "csfmask%mni#native.nii.gz"

                elif [  "${NormHandle}" == "n" ]; then
                  applywarp --ref="${T1_hr}.nii.gz" --in="${image_to_move}.nii.gz" \
                      --warp=t1hr2mni_NLwarp_inv --out="csfmask%mni#native.nii.gz" --interp=nn
                fi

                rm ${subj}_mask_csf_seg.nii.gz
                3dcalc -a ${T1_hr}_pve_0_thr${Seg_prob}.nii.gz -b "csfmask%mni#native.nii.gz" -expr "a*b" -prefix ${subj}_mask_csf_seg.nii.gz
              fi

              rm -f mask*.nii.gz

            done


            # checked already, via 3dROIstats -nzminmax -mask CoRe_011_mask_wm_invnorm_rs.nii.gz CoRe_011_mask_wm_invnorm_rs.nii.gz; wont interpolate numbers between 0 and 1 during resampling;

            rm ${subj}_mask_*_rs.nii.gz
            3dresample -master "${EPIbase}.nii.gz" -prefix "${subj}_mask_csf_${CSF_option}_rs.nii.gz" -inset "${subj}_mask_csf_${CSF_option}.nii.gz"
            3dresample -master "${EPIbase}.nii.gz" -prefix "${subj}_mask_wm_${WM_option}_rs.nii.gz" -inset "${subj}_mask_wm_${WM_option}.nii.gz"

            echo '|| '
            echo '|| ----- CHECK AGAIN, whether only 1 and 0 in the masks -----'
            echo '|| '
            echo "-------------------------------------------------------------------- "
            3dROIstats -nzminmax -mask "${subj}_mask_csf_${CSF_option}_rs.nii.gz" "${subj}_mask_csf_${CSF_option}_rs.nii.gz"
            echo "-------------------------------------------------------------------- "
            3dROIstats -nzminmax -mask "${subj}_mask_wm_${WM_option}_rs.nii.gz" "${subj}_mask_wm_${WM_option}_rs.nii.gz"
            echo "-------------------------------------------------------------------- "
            echo '|| '
            echo ' '


            # DO the coregistration step for sleep session
            if [ "${Pause_Dspike}" == "n" ]; then  # if Pause_Dspike=="y", chance are that "${data}Co.nii.gz" was already generated or faked inside sleep session, because the sleep session is so long,
                                                   # and sometime the reason to re-run this module is just for fine-tuning the normalization, so no need to repeat
              cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}
              rm "${data}Co.nii.gz"

              if [ "${Deoblique}" == "y" ]; then # coregistrate the sleep session to this plumb reference volume AGAIN.
                                                 # 2 do : no need to match center, because the base is from sleep session; WAIT, try match center first
                3dvolreg -base "${EPIbase}.nii.gz" -linear -twopass -noclip -coarserot -maxdisp -prefix "${data}Co.nii.gz" "${data}.nii.gz"
              else # pretend sleep session has been coregistrated within sleep session again; rename the sleep session for the simplicity;
                   # 2 do : whether should do twice?
                3dcalc -a "${data}.nii.gz" -expr "a" -prefix "${data}Co.nii.gz"
              fi
            fi


          # For other EPI sessions ;
        else  # simply just coregistrate other epi sessions with the sleep session base volume

            if [ "${Pause_Coreg}" == "n" ]; then # same logic as above
              cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}

              rm "${data}Co.nii.gz"
              # cp -f ${base_path}/${EPIbase}.nii.gz "${EPIbase}.nii.gz"
              # cp -f "${T1_path}/${T1_hr}.nii.gz" "${T1_hr}.nii.gz" do in next loop


              # match the centers first, if needed
              CMDist="$(@Center_Distance -dset ${data}.nii.gz ${EPIbase}.nii.gz)"
              echo
              echo "|| ----- CM distance = ${CMDist} mm-----"
              echo

              if [ `echo "${CMDist} > ${CMthreshold}" | bc` -eq 1 ]; then
                echo
                echo "|| ----- align the grid centers first, optionally possible between brain mass centers too -----"
                echo
                rm ${data}_bCMbackup.nii.gz
                3dcopy ${data}.nii.gz ${data}_bCMbackup.nii.gz
                @Align_Centers -base "${EPIbase}.nii.gz" -dset "${data}.nii.gz" -overwrite -no_cp #-prefix "${data}.nii.gz" #-cm #= brain's center ;By default, center = grid center.

              else
                echo
                echo "|| ----- =< ${CMthreshold} mm, no need to match center -----"
                echo
              fi



              # Do the Coregistration, EPI sessions to EPI reference image
              3dvolreg -base "${EPIbase}.nii.gz" -linear -twopass -noclip -coarserot -maxdisp -prefix "${data}Co.nii.gz" "${data}.nii.gz"

              # fine tuning Coregistration for some sbjs
              # 3dvolreg -base "${EPIbase}.nii.gz" -heptic -twopass -maxite 500 -x_thresh 0.003 -rot_thresh 0.007 -delta 1 -noclip -coarserot -maxdisp -prefix "${data}Co.nii.gz" "${data}.nii.gz"     # -edging 9         #
              # 3dvolreg -base "${EPIbase}.nii.gz" -linear -twopass -maxite 200 -x_thresh 0.003 -rot_thresh 0.007 -delta 0.2 -noclip -coarserot -maxdisp -verbose -verbose -prefix "${data}test.nii.gz" "${data}.nii.gz"
              # 3dvolreg -base "CoRe_220_coreg_base.nii" -linear -twopass -noclip -coarserot -maxdisp -prefix "test.nii.gz" "CoRe_220_D2_Retest_TSeq_MDs_bCMbackup.nii.gz"
            fi
          fi


        done

    # else
    fi

    #         ---- << code between loop >> ----

    ## [ SET UP ] : some vairables here, in case the above loops were skipped. just passing down some variables
    T1_path="${workdir}/${subj}/anat_HR"
    T1_hr="${subj}_anat_${T1_resolu}+BetCo"

    #         ---- << code between loop >> ----



    #   # ------------------------------------------------------------------------------------------------------------
    #   # ^w^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ <<  loop 3  >> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^o^
    if [ "${Loops2go}" == "all" ] || [ "${Loops2go}" == "3" ]; then

        for (( ss=0; ss<${sesslength}; ss++ ));      # because start from ss=0, and in bash index 0=1st item in the array
        do

          cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}
          data="${subj}_${daylist[$ss]}_${sessshort[$ss]}_MDsCo"


          echo ' '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo '|| '
          echo "||    -----[   4.  Noise removal: Trends(1/2) + Motion(6*4) + CSF/WM(2*4) = 33/34 parameters;  ]-----"
          echo "||    -----[   5.  Resampling with/without Normalization;  ]-----"
          echo "||    -----[   6.  Smoothing;  ]-----"
          echo "||    -----[   7.  Z-scoring;  ]-----"
          echo '|| '
          echo "||    current session under preprocessing is : [ '${sesslist[$ss]}' of ${subj} ]"
          echo "||                    the input file name is : [ ${data}.nii.gz ]"
          echo '|| '
          echo "-------------------------------------------------------------------------------------------------------------- "
          echo ' '


          # for checking convenience, copying the mni template into the folder
          rm -f *${T1_hr}* mnit1_09c_${T1_resolu}.nii.gz
          cp -f "${T1_path}/mnit1_09c_${T1_resolu}.nii.gz" "mnit1_09c_${T1_resolu}.nii.gz"


          # extract the averaged signal from physio-noise-masks
          rm ${subj}_*_mask_*.nii.gz
          cp -f "${T1_path}/${subj}_mask_csf_${CSF_option}_rs.nii.gz" "${subj}_mask_csf_${CSF_option}_rs.nii.gz"
          cp -f "${T1_path}/${subj}_mask_wm_${WM_option}_rs.nii.gz" "${subj}_mask_wm_${WM_option}_rs.nii.gz"

          for physio_mask in "${Masklist[@]}"
          do
            rm ${subj}_*_${physio_mask}.txt
            # rm ${subj}_${daylist[$ss]}_${sessshort[$ss]}_${physio_mask}.txt

            3dmaskave -quiet -mask "${subj}_mask_${physio_mask}_*_rs.nii.gz" \
            ${data}.nii.gz > ${subj}_${daylist[$ss]}_${sessshort[$ss]}_${physio_mask}.txt
          done


          ## < steps of doing nuissance regression >
            # 1. combine all regressors into one file; at this point, supposedly, we should have CSF, WM, and motion-6-columns;
            # 2. then the aim is to create their 1st order derivative(backward difference), quadratic term and combined with linear and quadratic(for sleep only) trends term.
                  # < ref paper > : https://doi.org/10.1016/j.neuroimage.2012.08.052
            # 3. do the regression


          # 1. combine all regressors into one file; at this point, supposedly, we should have CSF, WM, and motion-6-columns;
          rm all_regressors*.txt single_column_temp.txt
          1dcat ${subj}_*motion.txt ${subj}_*wm.txt ${subj}_*csf.txt > regressors_raw.txt

          unset rows scan
          scan=$(fslval ${data}.nii.gz dim4) # or afni version : scan=`3dinfo -nv ${data}.nii.gz`
          rows=$(more regressors_raw.txt|wc -l)


          if [ ${rows} -ne ${scan} ]; then
            echo ' '
            echo "|| ||||||||||||||||||||||||||||||||||||||||||||||||"
            echo "||   Time points do not match, or missing file(s) !"
            echo "||   Please check your data again! "
            echo "|| ||||||||||||||||||||||||||||||||||||||||||||||||"
            echo ' '

            exit 1
          fi

          unset cols nums
          nums=$(more regressors_raw.txt|wc -w)
          cols=$((nums / rows))



          # 2. then the aim is to create their 1st order derivative(backward difference), quadratic term and combined with linear and quadratic(for sleep only) trends term.
          1d_tool.py -infile "regressors_raw.txt" -set_nruns 1 -derivative -write "regressors_deriv.txt" # -derivative = -backward_diff
          1dcat regressors_raw.txt regressors_deriv.txt > all_regressors_temp.txt # combine them as temp for quadratic term calculation in next step
          1dcat regressors_raw.txt regressors_deriv.txt > regressors_raw+deriv.txt # used for Retest_TSeq, with only 34 time points

          # 1dmatcalc "&read(all_regressors_temp.txt) ^ &write(all_regressors.quadr.txt)"
          for (( cc=0; cc<${cols}*2; cc++ )); # quadratic term calculation and append to the end
          do
            1deval -a all_regressors_temp.txt[$cc] -expr 'a*a' > single_column_temp.txt
            1dcat all_regressors_temp.txt single_column_temp.txt > all_regressors_add1.txt

            rm all_regressors_temp.txt
            mv all_regressors_add1.txt all_regressors_temp.txt
            rm all_regressors_add1.txt

          done

          # the 3dTproject has the operation embedded; thus the code for creating two trends were commented
            # 1dcat `1deval -1D: -num ${scan} -expr 't'` `1deval -1D: -num ${scan} -expr 't^2'` > trends_2term.txt # create linear and quadratic trends;
            # 1dcat all_regressors_temp.txt trends_2term.txt > all_regressors.txt # combine all: 8+8+16+2= 34 parameters/regressors; 33 for normal sessions (not sleep no quadratic trend)

          # re-naming at the last
          mv all_regressors_temp.txt all_regressors.txt
          rm all_regressors_temp.txt single_column_temp.txt

          # 3. do the regression
          rm ${data}+Tavg.nii.gz BeforeNr_temp.nii.gz ${data}Nr.nii.gz
          fslmaths "${data}.nii.gz" -Tmean "${data}+Tavg.nii.gz"



          if [ ${scan} -ge 1000 ]   # as long the data is long enough > 1000;  if [ "${sessshort[$ss]}" == "sleep" ]; of notice there are multiple sleep sessions, even possible in the same day;
          then   # no need for the constant term; will be added automatically by 3dTproject; has checked via -verb operation
            3dTproject -automask -ort "all_regressors.txt" -prefix "BeforeNr_temp.nii.gz" -polort 2 -input "${data}.nii.gz" #-polort pp  = Remove polynomials up to and including degree pp. ++ Default value is 2; constant term will be added automatically
            echo
            echo "|| ----- ${data}.nii.gz has more than 1000 volumes, then quadratic trend also will be regressed out  -----"
            echo
          elif [ "${sessshort[$ss]}" == "Retest_TSeq" ]; then
            3dTproject -automask -ort "regressors_raw+deriv.txt" -prefix "BeforeNr_temp.nii.gz" -polort 1 -input "${data}.nii.gz" #-polort pp  = Remove polynomials up to and including degree pp. ++ Default value is 2; constant term will be added automatically
          else
            3dTproject -automask -ort "all_regressors.txt" -prefix "BeforeNr_temp.nii.gz" -polort 1 -input "${data}.nii.gz" #-polort pp  = Remove polynomials up to and including degree pp. ++ Default value is 2; constant term will be added automatically
          fi

          3dcalc -a "${data}+Tavg.nii.gz" -b "BeforeNr_temp.nii.gz" -expr "a+b" -prefix "${data}Nr.nii.gz"
          rm ${data}+Tavg.nii.gz BeforeNr_temp.nii.gz

          data="${data}Nr"

          ## 2do: check how prof.wu norm the anat and apply to epi
          echo ' '
          echo "-------------------------------------------------------------------- "
          echo '|| '
          echo "|| -----[   Resampling with/without Normalization;  ]-----"
          echo '|| '
          echo "-------------------------------------------------------------------- "
          echo ' '

          if [  "${Resample_Norm}" == "r" ]; then # only resampling no normalization

            rm ${data}Rs.nii.gz
            3dresample -dxyz ${Final_voxelsize} ${Final_voxelsize} ${Final_voxelsize} -prefix "${data}Rs.nii.gz" -input "${data}.nii.gz"
            data="${data}Rs"

          elif [ "${Resample_Norm}" == "rn" ]; then # normalization with resampling

            # prepare a final-voxel-size MNI template reference
            MNIref="${datapath}/Analysis_cabinet/mnit1_09c_1mm.nii.gz"

            rm MNIref.nii.gz
            if [  "${Final_voxelsize}" == "o" ]; then # keep it as original EPI voxel size
              3dZeropad -RL 100 -AP 100 -IS 80 -prefix ref4mni_ds.nii.gz CoRe_${subj}_coreg_base.nii.gz # dilate the box
              3dresample -master ref4mni_ds.nii.gz -prefix MNIref_temp.nii.gz -input ${MNIref}
              3dAutobox -input MNIref_temp.nii.gz -prefix MNIref_sbox.nii.gz -noclust # squeeze the box for later nonlinear normalization
              3dZeropad -RL 56 -AP 67 -IS 58 -prefix MNIref.nii.gz MNIref_sbox.nii.gz # recover the box. 3 numbers were estimated by prior try on direct downsampling. The way if with the -master is better, dont know why though.
              rm ref4mni_ds.nii.gz MNIref_*

              MNIref="MNIref.nii.gz"

            else # otherwise voxel size has been specified
              3dresample -dxyz ${Final_voxelsize} ${Final_voxelsize} ${Final_voxelsize} -prefix MNIref.nii.gz -input ${MNIref}
              MNIref="MNIref.nii.gz"

            fi


            # DO the normalization, the final voxel size will be inherited from the MNIref image
            if [  "${NormHandle}" == "l" ]; then # apply inv norm to physio_masks; either in a linear or nonlinear way
              rm xform_temp.mat xform_temp.nii.gz ${data}%mni.nii.gz
              flirt -ref "${MNIref}" -in "${data}.nii.gz" -applyxfm -init "${T1_path}/t1hr2mni.mat"\
                  -schedule ${FSLDIR}/etc/flirtsch/xyztrans.sch -usesqform \
                  -omat "xform_temp.mat" -out "xform_temp.nii.gz"
              flirt -ref "xform_temp.nii.gz" -in "${data}.nii.gz" -noresampblur -noresample -noclamp \
                  -applyxfm -init  "xform_temp.mat" -out "${data}%mni.nii.gz" # 2 do : check if resampled?

              data="${data}%mni"

            elif [  "${NormHandle}" == "n" ]; then
              rm ${data}%mniNL.nii.gz
              applywarp --ref="${MNIref}" --warp="${T1_path}/t1hr2mni_NLwarp.nii.gz" --in="${data}.nii.gz" --out="${data}%mniNL.nii.gz"

              # applywarp --ref="${datapath}/Analysis_cabinet/mnit1_09c_1mm.nii.gz" --in="CoRe_050_coreg_base.nii.gz" \
              #         --warp="${T1_path}/t1hr2mni_NLwarp.nii.gz" --out="CoRe_050_coreg_base%mniNL.nii.gz"

              data="${data}%mniNL"
            fi

          fi


          echo ' '
          echo "-------------------------------------------------------------------- "
          echo '|| '
          echo "|| -----[   spatial smoothing;  ]-----"
          echo '|| '
          echo "-------------------------------------------------------------------- "
          echo ' '

          # prepare a mask for below smoothing
          rm "Temp_anat_mask_temp.nii.gz"
          if [  "${Resample_Norm}" == "r" ] || [  "${Resample_Norm}" == "ne" ]; then # previously, with/without resampling, but no normalization; in native space, mask = T1 image after coregistration
            3dmask_tool -input "${T1_path}/${T1_hr}.nii.gz" -prefix "Temp_anat_mask_temp.nii.gz" -fill_holes -fill_dirs xy -dilate_input 1 -1
          elif [ "${Resample_Norm}" == "rn" ]; then # previously, done normalization with resampling, in MNI space, mask = MNI template
            3dmask_tool -input "${MNIref}" -prefix "Temp_anat_mask_temp.nii.gz" -fill_holes -fill_dirs xy -dilate_input 1 -1
          fi

          rm -f Temp_anat_mask.nii.gz ${data}S.nii.gz
          3dresample -master "${data}.nii.gz" -prefix "Temp_anat_mask.nii.gz" -inset "Temp_anat_mask_temp.nii.gz"
          rm "Temp_anat_mask_temp.nii.gz"

          # Do the smoothing
          3dBlurInMask -input "${data}.nii.gz" -FWHM "${Smth_kernel}" -mask Temp_anat_mask.nii.gz -float -quiet -prefix "${data}S.nii.gz"

          rm -f Temp_anat_mask.nii.gz

          ## using the automask to constrain the smoothing areas
          # 3dBlurToFWHM -input "${data}.nii.gz" -automask -FWHM "${Smth_kernel}" -detin -temper -prefix "${data}S.nii.gz"

          data="${data}S"

          echo ' '
          echo "-------------------------------------------------------------------- "
          echo '|| '
          echo "|| -----[   Z-scoring;  ]-----"
          echo '|| '
          echo "-------------------------------------------------------------------- "
          echo ' '

          # prepare the mean and std images
          rm -f ${data}+Tavg.nii.gz ${data}+Tstd.nii.gz ${data}Z.nii.gz

          3dTstat -mean -prefix "${data}+Tavg.nii.gz" "${data}.nii.gz"
          3dTstat -stdev -prefix "${data}+Tstd.nii.gz" "${data}.nii.gz"

          # z-score each time courses of the voxels
          3dcalc -a "${data}.nii.gz" -b "${data}+Tavg.nii.gz" -c "${data}+Tstd.nii.gz" -expr '(a-b)/c' -prefix "${data}Z.nii.gz"
          rm -f ${data}+Tavg.nii.gz ${data}+Tstd.nii.gz

          data="${data}Z"

          echo 'end 1 ' # for debugging convenience
      done
      echo 'end 2 '
    fi
    echo 'end 3 '




    # ------------------------------------------------------------------------------------------------------------
    # ^w^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ <<  loop checking  >> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^o^

    # copying all the resultant files into one directory for checking convenience
    if [ "${Loops2go}" == "all" ] || [ "${Loops2go}" == "chk" ]; then

      # create an independent checking-usage folder
      checkpath="MRI_PreP/CheckItOut"
      if [ ! -d "${datapath}/${checkpath}" ]; then
        mkdir "${datapath}/${checkpath}"
      fi
      cd "${datapath}/${checkpath}"
      if [ "${ReCheck}" == "y" ]; then
        rm -rf coreg2base/ norm2mni/ # saving more space by recreate the folder, checking subjects batch by batch
        mkdir "${datapath}/${checkpath}/coreg2base"
        mkdir "${datapath}/${checkpath}/norm2mni"
      fi

      # rm ${EPIbase}.nii.gz ${EPIbase}.nii.gz mnit1_09c_1mm.nii.gz

      ### copying target files into it
      rm ${datapath}/${checkpath}/coreg2base/${T1_hr}.nii.gz ${datapath}/${checkpath}/norm2mni/${T1_hr}%mni*.nii.gz
      3dcopy "${T1_path}/${T1_hr}.nii.gz" "${datapath}/${checkpath}/coreg2base/${T1_hr}.nii.gz" # also include the T1 coregistrated outcome
      if [  "${NormHandle}" == "l" ]; then
        3dcopy "${T1_path}/${T1_hr}%mni.nii.gz" "${datapath}/${checkpath}/norm2mni/${T1_hr}%mni.nii.gz" # also include the T1 normalization outcome
      elif [  "${NormHandle}" == "n" ]; then
        3dcopy "${T1_path}/${T1_hr}%mniNL.nii.gz" "${datapath}/${checkpath}/norm2mni/${T1_hr}%mniNL.nii.gz" # also include the T1 nonlinear normalization outcome
        3dcopy "${T1_path}/${T1_hr}%mni.nii.gz" "${datapath}/${checkpath}/norm2mni/${T1_hr}%mni.nii.gz" # also include linear normalization outcome for comparison

      fi


      for (( ss=0; ss<${sesslength}; ss++ ));      # because start from ss=0, and in bash index 0=1st item in the array
      do

        if [ "${ReCheck}" == "n" ]; then # just delete previous version if existed
          cd "${datapath}/${checkpath}"
          rm ${datapath}/${checkpath}/coreg2base/${subj}_${daylist[$ss]}*${sesslist[$ss]}*.nii.gz
          rm ${datapath}/${checkpath}/norm2mni/${subj}_${daylist[$ss]}*${sesslist[$ss]}*.nii.gz
        fi

        cd ${workdir}/${subj}/${daylist[$ss]}/${sesslist[$ss]}
        3dcopy `ls *Co.nii.gz` "${datapath}/${checkpath}/coreg2base/`ls *Co.nii.gz`"
        echo "[   > >> >>> >>>> >>>>> >>>>>>> >>>>>>>> >>>>>>>>> >>>>>>>>>> moving `ls *Co.nii.gz`   ]"

        if [  "${NormHandle}" == "l" ]; then
          3dcopy `ls *%mni.nii.gz` "${datapath}/${checkpath}/norm2mni/`ls *%mni.nii.gz`"
          echo "[   > >> >>> >>>> >>>>> >>>>>>> >>>>>>>> >>>>>>>>> >>>>>>>>>> moving `ls *%mni.nii.gz`   ]"
          # 3dcopy ${subj}_${daylist[$ss]}_${sessshort[$ss]}_MDsCoNr%mniSZ.nii.gz "${datapath}/${checkpath}/norm2mni/${subj}_${daylist[$ss]}_${sesslist[$ss]}_*%mni*Z.nii.gz"
          # echo "[   > >> >>> >>>> >>>>> >>>>>>> >>>>>>>> >>>>>>>>> >>>>>>>>>> moving ${subj}_${daylist[$ss]}_${sesslist[$ss]}_*%mni*Z.nii.gz   ]"

        elif [  "${NormHandle}" == "n" ]; then
          3dcopy `ls *%mniNL.nii.gz` "${datapath}/${checkpath}/norm2mni/`ls *MDsCoNr%mniNL.nii.gz`"
          echo "[   > >> >>> >>>> >>>>> >>>>>>> >>>>>>>> >>>>>>>>> >>>>>>>>>> moving `ls *%mniNL.nii.gz`   ]"
          # 3dcopy ${subj}_${daylist[$ss]}_${sessshort[$ss]}_MDsCoNr%mniNLSZ.nii.gz "${datapath}/${checkpath}/norm2mni/${subj}_${daylist[$ss]}_${sesslist[$ss]}_*%mniNL*Z.nii.gz"
          # echo "[   > >> >>> >>>> >>>>> >>>>>>> >>>>>>>> >>>>>>>>> >>>>>>>>>> moving ${subj}_${daylist[$ss]}_${sesslist[$ss]}_MDsCoNr%mniNLSZ.nii.gz   ]"

        fi

      done

      # additionally, collecting EPI base image for each subject
      3dcopy "${EPIbase}.nii.gz" "${datapath}/${checkpath}/coreg2base/${EPIbase}.nii.gz"

    fi
done

if [ "${Loops2go}" == "all" ] || [ "${Loops2go}" == "chk" ]; then
  # do not forget the latest reference images; take an advantage of the subjlist loop, to just do here just once
  3dcopy MNIref.nii.gz "${datapath}/${checkpath}/norm2mni/MNIref.nii.gz"
  3dcopy mnit1_09c_1mm.nii.gz "${datapath}/${checkpath}/norm2mni/mnit1_09c_1mm.nii.gz"
fi
