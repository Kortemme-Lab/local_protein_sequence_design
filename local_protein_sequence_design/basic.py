import pyrosetta
from pyrosetta import rosetta


def get_link_residues_task_op(sequence_symmetry_map):
    '''Generate the link residue task operation for a 
    given sequence_symmetry_map.
    '''
    # Initialize the equivalent residues

    equivalent_residues = []

    for source_seg in sequence_symmetry_map.keys():
        target_seg = sequence_symmetry_map[source_seg]

        for i in range(source_seg[1] - source_seg[0] + 1):
            equivalent_residues.append({source_seg[0] + i, target_seg[0] + i})

    # Merget the equivalent sets

    no_sets_to_merge = False

    while not no_sets_to_merge:
        no_sets_to_merge = True

        for i in range(len(equivalent_residues)): 
            for j in range(i + 1, len(equivalent_residues)):
                if len(equivalent_residues[i].intersection(equivalent_residues[j])) > 0:

                    equivalent_residues[i] = equivalent_residues[i].union(equivalent_residues[j])

                    del equivalent_residues[j]

                    no_sets_to_merge = False

                    break

            if not no_sets_to_merge:
                break

    # Create the link residues task operation

    lr = rosetta.protocols.task_operations.LinkResidues()
    
    for s in equivalent_residues:
        lr.add_group(','.join([str(r) for r in s]))
    
    return lr

def get_task_factory(pose, designable_residues, repackable_residues, extra_rotamers=True, limit_aro_chi2=True, layered_design=True, sequence_symmetry_map=None):
    '''Get a task factory given the designable and repackable residues.'''
    def list_to_str(l):
        return ','.join(list(str(i) for i in l))

    task_factory = rosetta.core.pack.task.TaskFactory()

    if len(designable_residues) > 0:
        designable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(designable_residues)) 
        racaa = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        racaa.aas_to_keep('GAPVILMFYWSTKRDENQ') # No CYS or HIS
        designable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                racaa, designable_selector)
        task_factory.push_back(designable_operation)

    if len(repackable_residues) > 0:
        repackable_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(repackable_residues)) 
        repackable_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.RestrictToRepackingRLT(), repackable_selector)
        task_factory.push_back(repackable_operation)

    natro_residues = [i for i in range(1, pose.size() + 1) if not i in designable_residues + repackable_residues]
    if len(natro_residues) > 0:
        natro_selector = rosetta.core.select.residue_selector.ResidueIndexSelector(list_to_str(natro_residues)) 
        natro_operation = rosetta.core.pack.task.operation.OperateOnResidueSubset(
                rosetta.core.pack.task.operation.PreventRepackingRLT(), natro_selector)
        task_factory.push_back(natro_operation)

    if extra_rotamers:
        ers = rosetta.core.pack.task.operation.ExtraRotamersGeneric()
        ers.ex1(True)
        ers.ex2(True)
        ers.extrachi_cutoff(18)

        task_factory.push_back(ers)

    if limit_aro_chi2:
        lac = rosetta.protocols.task_operations.LimitAromaChi2Operation()

        task_factory.push_back(lac)

    if layered_design:
        ld = rosetta.protocols.rosetta_scripts.XmlObjects.static_get_task_operation(
            '''<LayerDesign name="layer_all" layer="core_boundary_surface_Nterm_Cterm" use_sidechain_neighbors="True">
    		<Nterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Nterm>
    		<Cterm>
    			<all append="DEGHKNQRST" />
    			<all exclude="CAFILMPVWY" />
    		</Cterm>
        </LayerDesign>''')
        task_factory.push_back(ld)

#    if sequence_symmetry_map:
#        task_factory.push_back(get_link_residues_task_op(sequence_symmetry_map))

    return task_factory

def set_non_crystal_symmetry(pose, score_function, sequence_symmetry_map):
    '''Set up non crystal symmetry for a pose, a score function
    given a sequence_symmetry_map. 
    A sequence_symmetry_map is defined as a dictionary 
    {(target_start, target_stop) : (source_start, source_stop)}
    '''
    # Set up the constraints
    
    ncs_mover = rosetta.protocols.symmetry.SetupNCSMover()

    ncs_mover.set_bb(True)
    ncs_mover.set_chi(False)
    ncs_mover.set_weight(0.01)
    ncs_mover.set_limit(10)

    for target in sequence_symmetry_map.keys():
        source = sequence_symmetry_map[target]
        ncs_mover.add_group('{0}-{1}'.format(source[0], source[1]), '{0}-{1}'.format(target[0], target[1]))

    ncs_mover.apply(pose)

    # Update the score function weights

    score_function.set_weight(rosetta.core.scoring.dihedral_constraint, 10)
    score_function.set_weight(rosetta.core.scoring.atom_pair_constraint, 1)
