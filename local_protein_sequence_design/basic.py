import pyrosetta
from pyrosetta import rosetta


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

    if sequence_symmetry_map:
        task_factory.push_back(get_link_residues_task_operation(sequence_symmetry_map, pose.size()))

    return task_factory

def get_link_residues_task_operation(sequence_symmetry_map, pose_size):
    '''Get a LinkResidues task operation for a sequence_symmetry_map.
    A sequence_symmetry_map is defined as a dictionary 
    {(target_start, target_stop) : (source_start, source_stop)}
    
    Note: Although I implementated this function. It seems not useful
    because LinkResidues task operation usually causes seg-fault.
    I think that the problem is due to two linked residues have different
    allowed residues.
    '''
    # Get the equivalent groups.
    # Note the current implementation only works if no source residue
    # is also a target residue 

    equivalent_groups = {}

    for target_seg in sequence_symmetry_map.keys():
        target_start = target_seg[0]
        target_stop = target_seg[1]
        source_start = sequence_symmetry_map[target_seg][0]
        source_stop = sequence_symmetry_map[target_seg][1]

        for i in range(target_stop - target_start + 1):
            t = target_start + i
            s = source_start + i

            if s in equivalent_groups.keys():
                equivalent_groups[s].append(t)
            else:
                equivalent_groups[s] = [s, t]

    # Create the LinkResidues task operation

    link_residues = rosetta.protocols.task_operations.LinkResidues()

    for k in equivalent_groups.keys():
        link_residues.add_group(','.join(str(x) for x in equivalent_groups[k] if x != 1 and x != pose_size))

    link_residues.add_group('1,1')
    link_residues.add_group('{0},{0}'.format(pose_size))

    return link_residues
