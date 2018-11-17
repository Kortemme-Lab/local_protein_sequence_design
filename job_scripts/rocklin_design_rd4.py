#!/usr/bin/env python3
"""run design using the rosetta scripts protocol for round 4 of EHEE topology in
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
        <Add filter_name="cavity_threshold" /> ?
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

xmlobj = rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
'''
    <SCOREFXNS>
        <ScoreFunction name="sfxn_std"  weights="ref2015">
            <Reweight scoretype="aa_composition" weight="1.0" />
            <Reweight scoretype="p_aa_pp" weight="0.6" /> #from 0.4 
        </ScoreFunction>
        <ScoreFunction name="nov15"  weights="ref2015"/>
        <ScoreFunction name="net_atr_net_sol"  weights="empty.wts">
            <Reweight scoretype="fa_atr" weight="1.0" />
            <Reweight scoretype="fa_rep" weight="0.55" /> 
            <Reweight scoretype="fa_sol" weight="1.0" />
            <Reweight scoretype="fa_elec" weight="1.0" />
        </ScoreFunction>
        <ScoreFunction name="SFXN1" weights="fldsgn_cen">
            <Reweight scoretype="hbond_sr_bb" weight="1.0" />
            <Reweight scoretype="hbond_lr_bb" weight="1.0" />
            <Reweight scoretype="atom_pair_constraint" weight="1.0" />
            <Reweight scoretype="angle_constraint" weight="1.0" />
            <Reweight scoretype="dihedral_constraint" weight="1.0" />
        </ScoreFunction>
        <ScoreFunction name="SFXN2" weights="fldsgn_cen">
            <Reweight scoretype="hbond_sr_bb" weight="1.0" />
            <Reweight scoretype="hbond_lr_bb" weight="1.0" />
            <Reweight scoretype="atom_pair_constraint" weight="1.0" />
            <Reweight scoretype="angle_constraint" weight="1.0" />
            <Reweight scoretype="dihedral_constraint" weight="1.0" />
        </ScoreFunction>
        <ScoreFunction name="TotalHydrophobic" weights="total_hydrophobic_weights.wts"/>
        <ScoreFunction name="reg_nov15"  weights="ref2015">
            <Reweight scoretype="aa_composition" weight="1.0" />
            <Reweight scoretype="p_aa_pp" weight="0.7" /> #from 0.6
        </ScoreFunction>
    </SCOREFXNS>
    
    <TASKOPERATIONS>
        <LayerDesign name="layer_core_SCN" layer="core" pore_radius="2.0" verbose="true" use_sidechain_neighbors="True" core="4" />
        <ReadResfile name="resfile" filename="ehee.resfile"/>
    </TASKOPERATIONS>
    
    <RESIDUE_SELECTORS>
        <ResidueName name="big_np" residue_name3="ALA,CYS,ASP,GLU,GLY,HIS,LYS,MET,ASN,PRO,GLN,ARG,SER,THR,VAL"/>
    </RESIDUE_SELECTORS>
    
    <FILTERS>
        <ResidueCount name="res_count_all" max_residue_count="9999" confidence="0"/>
        <ResidueCount name="np_count" residue_types="PHE,MET,ILE,LEU,TYR,TRP,VAL" min_residue_count="9" confidence="1"/>
        <ScoreType name="hbond_sfn" scorefxn="sfxn_std" score_type="hbond_lr_bb" threshold="0"/>
        <CalculatorFilter name="bb" equation="hbond / rescount" threshold="-0.30" confidence="1">
            <Var name="hbond" filter="hbond_sfn"/>
            <Var name="rescount" filter="res_count_all"/>
        </CalculatorFilter>
        <HelixKink name="hk1" blueprint="ehee.bp"/>
        <SheetTopology name="sf1" blueprint="ehee.bp" />
        <SecondaryStructure name="ss1" blueprint="ehee.bp.ss" />
        <CompoundStatement name="cs1">
            <AND filter_name="ss1" />
            <AND filter_name="hk1" />
            <AND filter_name="sf1" />
        </CompoundStatement>
        <ScoreType name="total_score_cen" score_type="total_score" scorefxn="SFXN2" confidence="0" threshold="0" />
        <ScoreType name="net_atr_net_sol" score_type="total_score" scorefxn="net_atr_net_sol" confidence="0" threshold="0" />
        <ScoreType name="omega" score_type="omega" scorefxn="sfxn_std" confidence="1" threshold="2.5" />
        <ScoreType name="p_aa_pp" score_type="p_aa_pp" scorefxn="sfxn_std" confidence="1" threshold="-12" />
        <AverageDegree name="degree" confidence="1" threshold="9.2"/>
        <ExposedHydrophobics name="exposed" confidence="0"/>
        <AtomicContactCount name="contact" confidence="1"/>
        <CalculatorFilter name="contact_threshold" equation="-c" threshold="-125" confidence="1" >
            <Var name="c" filter="contact"/>
        </CalculatorFilter>

        <ResidueCount name="AlaCount" residue_types="ALA" max_residue_count="6" confidence="1"/>
        <AverageDegree name="degree_core_SCN" task_operations="layer_core_SCN" confidence="1" threshold="9.4" />
        <ResidueCount name="res_count_core_SCN" task_operations="layer_core_SCN" max_residue_count="9999" confidence="0"/>
        <CalculatorFilter name="percent_core_SCN" equation="- rescount_coreSCN / (rescount3 + 0.01)" threshold="-0.1" confidence="1" >
            <Var name="rescount3" filter="res_count_all"/>
            <Var name="rescount_coreSCN" filter="res_count_core_SCN"/>
        </CalculatorFilter>
        <BuriedUnsatHbonds2 name="unsat_hbond" confidence="1" jump_number="0" cutoff="1"/>
        <TotalSasa name="exposed_hydrophobics" confidence="1" hydrophobic="True" upper_threshold="1700"/>
        <ScoreType name="total_hydrophobic" scorefxn="TotalHydrophobic" threshold="0" />
        <CalculatorFilter name="buried_np" equation="total - exposed" threshold="1" confidence="0">
            <Var name="total" filter="total_hydrophobic"/>
            <Var name="exposed" filter="exposed_hydrophobics"/>
        </CalculatorFilter>
        <CalculatorFilter name="buried_np_per_res" equation="-buried / res" threshold="-118" confidence="1">
            <Var name="buried" filter="buried_np"/>
            <Var name="res" filter="res_count_all"/>
        </CalculatorFilter>
        <ScoreType name="total_score" scorefxn="nov15" threshold="0"/>
        <ScoreType name="fa_atr_filter" scorefxn="nov15" threshold="0" score_type="fa_atr" />
        <CalculatorFilter name="score_per_res" equation="total_score / res" threshold="-3.2" confidence="1">
            <Var name="total_score" filter="total_score"/>
            <Var name="res" filter="res_count_all"/>
        </CalculatorFilter>
        <CalculatorFilter name="fa_atr_per_res" equation="fa_atr_score / res" threshold="-5.2" confidence="1">
            <Var name="fa_atr_score" filter="fa_atr_filter"/>
            <Var name="res" filter="res_count_all"/>
        </CalculatorFilter>
        <CalculatorFilter name="net_atr_net_sol_per_res" equation="net_atr_net_sol_score / res" threshold="-3.15" confidence="1">
            <Var name="net_atr_net_sol_score" filter="net_atr_net_sol"/>
            <Var name="res" filter="res_count_all"/>
        </CalculatorFilter>
        <SSPrediction name="mismatch_probability" confidence="1" cmd="runpsipred_single" use_probability="1" mismatch_probability="True" use_svm="0" threshold="0.34"/>
    </FILTERS>
    
    <TASKOPERATIONS>
        <LimitAromaChi2 name="limitchi2" include_trp="1" />
        <LayerDesign name="layer_surfacedes" layer="surface" surface_E="90" surface_H="60" pore_radius="1.8">
            <surface>
                <all append="ADEHIKQRTV"/>
                <all exclude="CFGLMNPSWY" />
                <Helix append="ADEHKQR" />
                <Helix exclude="CFGILMNPSTVWY" />
                <Strand append="ADEHIKNQRTV" />
                <Strand exclude="CFGLMPSWY" />
            </surface>
        </LayerDesign>
        <OperateOnResidueSubset name="BigNPOnly" selector="big_np">
            <RestrictToRepackingRLT/>
        </OperateOnResidueSubset>
    </TASKOPERATIONS>
    
    <MOVERS>
        <Dssp name="dssp" />
        <SheetCstGenerator name="sheet_new1" cacb_dihedral_tolerance="0.6" blueprint="ehee.bp" />
        <SetSecStructEnergies name="set_ssene1" scorefxn="SFXN1" blueprint="ehee.bp" />	
        <BluePrintBDR name="bdr1" use_abego_bias="1" scorefxn="SFXN1" constraint_generators="sheet_new1" constraints_NtoC="0" blueprint="ehee.bp" />
        <DumpPdb name="dump" fname="pass" tag_time="True"/>
    <DumpPdb name="dump2" fname="pre_sf" tag_time="True"/>
        <FastDesign name="fastdes" task_operations="limitchi2,resfile" scorefxn="sfxn_std" clear_designable_residues="0" repeats="4" ramp_down_constraints="0" />
        <FastDesign name="fastdes_surface" task_operations="limitchi2,layer_surfacedes,BigNPOnly,resfile" scorefxn="sfxn_std" clear_designable_residues="0" repeats="1" ramp_down_constraints="0" />
        <ParsedProtocol name="build_dssp1" >
            <Add mover_name="bdr1" />
            <Add mover_name="dssp" />
            <Add filter_name="ss1" />
            <Add filter_name="hk1" />
            <Add filter_name="degree" />
            <Add mover_name="fastdes"/>
            <Add mover_name="fastdes_surface" />
            <Add filter_name="np_count" />   
            <Add filter_name="AlaCount" />   
            <Add filter_name="unsat_hbond" /> #?
            <Add filter_name="buried_np_per_res" /> 
            <Add filter_name="score_per_res" /> 
            <Add filter_name="net_atr_net_sol_per_res" />
            <Add filter_name="exposed_hydrophobics" /> 
            <Add filter_name="mismatch_probability" />
            <Add mover_name="dump"/>
        </ParsedProtocol>
        <LoopOver name="lover1" mover_name="build_dssp1" iterations="500" drift="0" ms_whenfail="FAIL_DO_NOT_RETRY" />
        <ParsedProtocol name="phase1" >
            <Add mover_name="set_ssene1" />
            <Add mover_name="lover1" />
        </ParsedProtocol>
    </MOVERS>
''')

# PROTOCOLS
# <Add mover_name="phase1" />
# <Add mover_name="dssp" />
# <Add filter_name="total_score_cen" />
# <Add filter_name="cs1" />
phase1 = xmlobj.get_mover('phase1')
dssp = xmlobj.get_mover('dssp')
total_score_cen = xmlobj.get_filter('total_score_cen')
cs1 = xmlobj.get_filter('cs1')

pose = pose_from_file('start.pdb')

phase1.apply(pose)
dssp.apply(pose)
total_score_cen.apply(pose)
cs1.apply(pose)
# pose.dump_pdb('design.pdb.gz')
