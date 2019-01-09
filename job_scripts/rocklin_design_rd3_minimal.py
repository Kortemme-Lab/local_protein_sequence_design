#!/usr/bin/env python3
"""run design using the rosetta scripts protocol for round 3 of EHEE topology in
Rocklin, Gabriel J., et al.
"Global analysis of protein folding using massively parallel design, synthesis, and testing."
Science 357.6347 (2017): 168-175.

modifications made in order to allow design for EHEE topologies of different length:
replaced all beta_* score fxns with current default - ref2015
all uses of CavityVolume - causes an error
    <CavityVolume name="cavity" confidence="0"/>
    <CalculatorFilter name="cavity_threshold" equation="c" threshold="25" confidence="1" >
        <Var name="c" filter="cavity"/>
    </CalculatorFilter>
    <Add filter_name="cavity_threshold" />

made blueprint movers and filters use design.bp which is made from dssp of input structure
replaced BuriedUnsatHbonds2 with BuriedUnsatHbonds
"""

from pyrosetta import *

init(options='-aa_composition_setup_file ehee.hydrophobic.comp '
             '-ex1 -ex2aro -use_input_sc -no_his_his_pairE -nblist_autoupdate true '
             '-chemical:exclude_patches LowerDNA UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB '
             'VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 '
             'ser_phosphorylated thr_phosphorylated tyr_phosphorylated tyr_sulfated '
             'lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated '
             'cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm '
             '-relax::ramp_constraints true '
             '-mute core.pack.interaction_graph.interaction_graph_factory'
             '-mute core.scoring.rms_util '
             '-mute core.pack.task '
             '-mute core.scoring.NeighborList '
             '-mute core.pack.annealer.MultiCoolAnnealer '
             '-mute core.pack.pack_rotamers '
             '-mute protocols.forge.remodel.RemodelDesignMover')

pdb_name = 'EHEE_rd1_0284.pdb'
pose = pose_from_file(pdb_name)
dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()
seq_str = pose.sequence()

blueprint_lines = ['SSPAIR 1-3.A.0;2-3.A.0\n']
for index, dssp_char in enumerate(dssp_str):
    amino_acid = seq_str[index]
    if dssp_char == 'E':
        blueprint_ss = 'EB'
    elif dssp_char == 'H':
        blueprint_ss = 'HA'
    elif dssp_char == 'L':
        blueprint_ss = 'LX'
    else:
        raise Exception('Unexpected dssp character {0} in {1}'.format(dssp_char, pdb_name))
    blueprint_lines.append('{0}    {1}    {2}    .\n'.format(index + 1, amino_acid, blueprint_ss))

with open('design.bp', 'w') as o:
    o.writelines(blueprint_lines)

xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''
    <SCOREFXNS>
        <ScoreFunction name="sfxn_std"  weights="ref2015">
          <Reweight scoretype="aa_composition" weight="1.0" />
          <Reweight scoretype="p_aa_pp" weight="0.8" /> #from 0.4 
          <Reweight scoretype="rama" weight="0.30" />  #from 0.25
          <Reweight scoretype="fa_elec" weight="0.5" /> #from 1.0
    
        </ScoreFunction>
        <ScoreFunction name="ref2015"  weights="ref2015"/>
        <ScoreFunction name="SFXN1" weights="fldsgn_cen">
            #Reweight scoretype="cenpack" weight="1.0" />
            <Reweight scoretype="hbond_sr_bb" weight="1.0" />
            <Reweight scoretype="hbond_lr_bb" weight="1.0" />
            <Reweight scoretype="atom_pair_constraint" weight="1.0" />
            <Reweight scoretype="angle_constraint" weight="1.0" />
            <Reweight scoretype="dihedral_constraint" weight="1.0" />
        </ScoreFunction>
    
        <ScoreFunction name="SFXN2" weights="fldsgn_cen">
            #Reweight scoretype="cenpack" weight="1.0" />
            <Reweight scoretype="hbond_sr_bb" weight="1.0" />
            <Reweight scoretype="hbond_lr_bb" weight="1.0" />
            <Reweight scoretype="atom_pair_constraint" weight="1.0" />
            <Reweight scoretype="angle_constraint" weight="1.0" />
            <Reweight scoretype="dihedral_constraint" weight="1.0" />
        </ScoreFunction>
    
        <ScoreFunction name="TotalHydrophobic" weights="total_hydrophobic_weights.wts"/>
    
    </SCOREFXNS>
    <TASKOPERATIONS>
         <LayerDesign name="layer_core_SCN" layer="core" pore_radius="2.0" verbose="true" use_sidechain_neighbors="True" core="4" make_pymol_script="1"/>
    </TASKOPERATIONS>
    
    
    <FILTERS>
          <ResidueCount name="res_count_all" max_residue_count="9999" confidence="0"/>
        <ResidueCount name="np_count" residue_types="PHE,MET,ILE,LEU,TYR,TRP,VAL" min_residue_count="7" confidence="0"/>
        <ResidueCount name="np_count_ala" residue_types="PHE,MET,ILE,LEU,TYR,TRP,VAL,ALA" min_residue_count="11" confidence="0"/>
        <ScoreType name="hbond_sfn" scorefxn="sfxn_std" score_type="hbond_lr_bb" threshold="0"/>
    
        <CalculatorFilter name="bb" equation="hbond / rescount" threshold="-0.2" confidence="0">
            <Var name="hbond" filter="hbond_sfn"/>
            <Var name="rescount" filter="res_count_all"/>
        </CalculatorFilter>
    
    
        <HelixKink name="hk1" blueprint="design.bp"/>
        <SheetTopology name="sf1" blueprint="design.bp"/>
        <SecondaryStructure name="ss1" blueprint="design.bp"/> 
        <CompoundStatement name="cs1">
            <AND filter_name="ss1" />
            <AND filter_name="hk1" />
            <AND filter_name="sf1"/> 
        </CompoundStatement>
    
    
        <ScoreType name="total_score_cen" score_type="total_score" scorefxn="SFXN2" confidence="0" threshold="0" />
    
        <AverageDegree name="degree" confidence="1" threshold="9.25"/>
        <ExposedHydrophobics name="exposed" confidence="0"/>
        <AtomicContactCount name="contact" confidence="1"/>
    
        <CalculatorFilter name="contact_threshold" equation="-c" threshold="240" confidence="0" >
                <Var name="c" filter="contact"/>
            </CalculatorFilter>
    
    
            <ResidueCount name="AlaCount" residue_types="ALA" max_residue_count="8" confidence="0"/>
            <AverageDegree name="degree_core_SCN" task_operations="layer_core_SCN" confidence="0" threshold="9.4" />
    
    
            <ResidueCount name="res_count_core_SCN" task_operations="layer_core_SCN" max_residue_count="9999" confidence="0"/>
    
        <CalculatorFilter name="percent_core_SCN" equation="- rescount_coreSCN / (rescount3 + 0.01)" threshold="-0.1" confidence="0" >
                <Var name="rescount3" filter="res_count_all"/>
                <Var name="rescount_coreSCN" filter="res_count_core_SCN"/>
            </CalculatorFilter>
    
            <BuriedUnsatHbonds name="unsat_hbond" confidence="0" jump_number="0" cutoff="5"/>
    
          <TotalSasa name="exposed_hydrophobics" confidence="0" hydrophobic="True" upper_threshold="1650"/>
          <ScoreType name="total_hydrophobic" scorefxn="TotalHydrophobic" threshold="0" />
    
          <CalculatorFilter name="buried_np" equation="total - exposed" threshold="1" confidence="0">
              <Var name="total" filter="total_hydrophobic"/>
              <Var name="exposed" filter="exposed_hydrophobics"/>
          </CalculatorFilter>
    
          <CalculatorFilter name="buried_np_per_res" equation="-buried / res" threshold="-116" confidence="0">
              <Var name="buried" filter="buried_np"/>
              <Var name="res" filter="res_count_all"/>
          </CalculatorFilter>
    
        <ScoreType name="total_score" scorefxn="ref2015" threshold="0"/>
        <ScoreType name="fa_atr_filter" scorefxn="ref2015" threshold="0" score_type="fa_atr" />
    
    
          <CalculatorFilter name="score_per_res" equation="total_score / res" threshold="-2.1" confidence="0">
              <Var name="total_score" filter="total_score"/>
              <Var name="res" filter="res_count_all"/>
          </CalculatorFilter>
    
          <CalculatorFilter name="fa_atr_per_res" equation="fa_atr_score / res" threshold="-4.125" confidence="0">
              <Var name="fa_atr_score" filter="fa_atr_filter"/>
              <Var name="res" filter="res_count_all"/>
          </CalculatorFilter>
    
      <PackStat name="pack" confidence="0" threshold="0.6"/>
      SSShapeComplementarity name="ss_sc" verbose="0" confidence="0" min_sc="0.67"/>
        <SSPrediction name="mismatch_probability" confidence="0" cmd="runpsipred_single" use_probability="1" mismatch_probability="True" use_svm="0" threshold="0.35"/>
    </FILTERS>
    <TASKOPERATIONS>
        <LimitAromaChi2 name="limitchi2" include_trp="1" />
        <LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True" pore_radius="2.0" verbose="true" core="3.5" surface="1.4" make_pymol_script="1" >
               <Nterm>
               <all append="DEGHKNQRST" />
               <all exclude="CAFILMPVWY" />
             </Nterm>
             <Cterm>
               <all append="DEGHKNQRST" />
               <all exclude="CAFILMPVWY" />
            </Cterm>
          </LayerDesign>
    
    
        <NoRepackDisulfides name="exemptdisulf" />
    </TASKOPERATIONS>
    <MOVERS>
        <Dssp name="dssp" />
    
        <SheetCstGenerator name="sheet_new1" cacb_dihedral_tolerance="0.6" blueprint="design.bp"/> 
        <SetSecStructEnergies name="set_ssene1" scorefxn="SFXN1" blueprint="design.bp"  />
        <BluePrintBDR name="bdr1" use_abego_bias="1" scorefxn="SFXN1" constraint_generators="sheet_new1" constraints_NtoC="0" blueprint="design.bp"/>
        <DumpPdb name="dump" fname="pass" tag_time="True"/>
    
    
        <FastDesign name="fastdes" task_operations="limitchi2,layer_all" scorefxn="sfxn_std" clear_designable_residues="0" repeats="7ƒb" ramp_down_constraints="0" />
        <FastDesign name="fastdes4" task_operations="limitchi2,layer_all" scorefxn="sfxn_std" clear_designable_residues="0" repeats="4" ramp_down_constraints="0" />
    
    
        <ParsedProtocol name="build_dssp1" >
            #Add mover_name="bdr1" />
            <Add mover_name="dssp"/> 
            <Add filter_name="cs1"/>
            <Add filter_name="degree"/> 
            <Add mover_name="fastdes"/>
            <Add filter_name="np_count_ala" />
            <Add filter_name="contact_threshold" />
            <Add filter_name="unsat_hbond" />
            <Add filter_name="pack" />
            <Add filter_name="buried_np_per_res" />
            <Add filter_name="score_per_res" />
            <Add filter_name="fa_atr_per_res" />
            <Add filter_name="exposed_hydrophobics" />
            #Add filter_name="mismatch_probability" />
            <Add mover_name="dump"/>
    
        </ParsedProtocol>
        <LoopOver name="lover1" mover_name="build_dssp1" iterations="1" drift="0" ms_whenfail="FAIL_DO_NOT_RETRY" />
        <ParsedProtocol name="phase1" >
            <Add mover_name="set_ssene1" />
            <Add mover_name="lover1" />
        </ParsedProtocol>
    
    </MOVERS>
    ''')

# PROTOCOLS
# <Add mover_name="phase1" />
phase1 = xmlobj.get_mover('phase1')
phase1.apply(pose)

# fastdes = xmlobj.get_mover('fastdes')
# fastdes.apply(pose)
# pose.dump_pdb('test_minimal.pdb')
