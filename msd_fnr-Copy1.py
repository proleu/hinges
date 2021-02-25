import pyrosetta
import os
import glob
import shutil
import sys

index = int(sys.argv[1])
np_penalty =str(sys.argv[2])

pyrosetta.init("-dunbrack_prob_buried 0.9 -dunbrack_prob_nonburied 0.9 -dunbrack_prob_buried_semi 0.9 -dunbrack_prob_nonburied_semi 0.9 -beta_nov16")
sf_obj = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
    '''<SCOREFXNS>
        <ScoreFunction name="sfxn" weights="beta_nov16" symmetric="1" />
        <ScoreFunction name="sfxn_design" weights="beta_nov16" symmetric="1" >
            <Reweight scoretype="res_type_constraint" weight="2.0" />
            <Reweight scoretype="aa_composition" weight="1.0" />
            <Set use_hb_env_dep="true" />
            #Reweight scoretype="hbnet" weight="0.7" />
            #Set elec_context_dependent="true" /> # not in master yet
            <Reweight scoretype="approximate_buried_unsat_penalty" weight="17" />
            <Set approximate_buried_unsat_penalty_burial_atomic_depth="3.5" />
            <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-1.0" />
            <Set approximate_buried_unsat_penalty_natural_corrections1="true" />
            <Set approximate_buried_unsat_penalty_hbond_bonus_cross_chain="-7" />
            <Set approximate_buried_unsat_penalty_hbond_bonus_ser_to_helix_bb="1"/>
            # lk_ball is slooooooooooow
            <Reweight scoretype="lk_ball" weight="0" />
            <Reweight scoretype="lk_ball_iso" weight="0" />
            <Reweight scoretype="lk_ball_bridge" weight="0" />
            <Reweight scoretype="lk_ball_bridge_uncpl" weight="0" />
        </ScoreFunction>

    </SCOREFXNS>''')

sf = sf_obj.get_score_function('sfxn_design')
sf_clean = sf_obj.get_score_function('sfxn')


file_x = glob.glob("0*pdb")[0]
file_y = glob.glob("1*pdb")[0]


pose1 = pyrosetta.io.pose_from_pdb(file_x)
pose2 = pyrosetta.io.pose_from_pdb(file_y)
#setup_symm = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover("C2_Z.sym")
#setup_symm.apply(pose1)
#setup_symm.apply(pose2)
num_des_res = []
des_res = []
score_list = []
ala_penalty = 1


def des_y_initial(despose,refpose,weight=0,strict_layers = 0):
    true_sel = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    allres = pyrosetta.rosetta.core.select.get_residues_from_subset(true_sel.apply(despose))
    diff = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for i in allres:
        if despose.sequence(i,i) == "C": # maintain disulfides in despose
            continue
        elif refpose.sequence(i,i) == "C": # safely replace despose residue with CYS (not CYD)
            mut = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
            mut.set_target(i)
            mut.set_res_name(pyrosetta.rosetta.core.chemical.AA(2)) # 2 is CYS
            mut.apply(despose)
        elif despose.sequence(i,i) != refpose.sequence(i,i):
            diff.append(i)
            despose.replace_residue(i,refpose.residue(i),1)
        else:
            pass
    num_des_res.append(len(diff))
    des_res.append(diff)
    #what to design
    objs_sse = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        '''
        <RESIDUE_SELECTORS>
            <SSElement name="part1" selection="n_term" to_selection="4,H,E" chain="A" reassign_short_terminal_loop="2" />
            <SSElement name="part2" selection="5,H,S" to_selection="c_term" chain="A" reassign_short_terminal_loop="2" />
        </RESIDUE_SELECTORS>
        ''')
    part1 = objs_sse.get_residue_selector('part1')
    part2 = objs_sse.get_residue_selector('part2')
    intsel = pyrosetta.rosetta.core.select.residue_selector.InterGroupInterfaceByVectorSelector(part1,part2)
    intdes = pyrosetta.rosetta.core.select.get_residues_from_subset(intsel.apply(despose))
    intref = pyrosetta.rosetta.core.select.get_residues_from_subset(intsel.apply(refpose))
    intall = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for i in intdes:
        intall.append(i)
    for i in intref:
        intall.append(i)
    designable = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(intall)
    #designable = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(diff)
    packable = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(designable,6,True)
    pack_option = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
    pack = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pack_option, designable, True)
    lock_option = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    lock = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(lock_option, packable, True)
    #some task ops
    aroChi = pyrosetta.rosetta.protocols.task_operations.LimitAromaChi2Operation()
    aroChi.chi2max(110)
    aroChi.chi2min(70)
    aroChi.include_trp(True)
    #layer design - selectors -- 
    ss1 = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose1)
    ss2 = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose2)
    surf_sel = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    surf_sel.set_layers(0,0,1)
    surf_sel.set_use_sc_neighbors(0)
    surf_sel.set_cutoffs(20,50)
    surf1 = pyrosetta.rosetta.core.select.get_residues_from_subset(surf_sel.apply(pose1))
    surf2 = pyrosetta.rosetta.core.select.get_residues_from_subset(surf_sel.apply(pose2))
    core_sel = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    core_sel.set_layers(1,0,0)
    core_sel.set_use_sc_neighbors(0)
    core1 = pyrosetta.rosetta.core.select.get_residues_from_subset(core_sel.apply(pose1))
    core2 = pyrosetta.rosetta.core.select.get_residues_from_subset(core_sel.apply(pose2))
    core_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    surf_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_core = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_surf = pyrosetta.rosetta.utility.vector1_unsigned_long()
    surf_core = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for i in allres:
        if i in core1:
            if i in core2:
                core_both.append(i)
            elif i in surf2:
                surf_core.append(i)
            else:
                bdry_core.append(i)
        elif i in surf1:
            if i in surf2:
                surf_both.append(i)
            elif i in core2:
                surf_core.append(i)
            else:
                bdry_surf.append(i)
        else:
            if i in core2:
                bdry_core.append(i)
            elif i in surf2:
                bdry_surf.append(i)
            else:
                bdry_both.append(i)
    if len(core_both) > 0:
        sel_core_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(core_both)
    else:
        sel_core_both = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    sel_surf_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(surf_both)
    if len(bdry_core) > 0:
        sel_bdry_core = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_core)
    else:
        sel_bdry_core = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    if len(bdry_surf) > 0:
        sel_bdry_surf = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_surf)
    else:
        sel_bdry_surf = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    if len(surf_core) > 0:
        sel_surf_core = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(surf_core)
    else:
        sel_surf_core = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    sel_bdry_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_both)
    if strict_layers == 1:
        sel_c = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_core_both,sel_bdry_core)
        sel_b = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_bdry_both,sel_surf_core)
        sel_s = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_surf_both,sel_bdry_surf)
    else:
        sel_c = sel_core_both
        sel_s = sel_surf_both
        sel_c_or_s = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_core_both,sel_surf_both)
        sel_b = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(sel_c_or_s)
        
    objs_sel = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        '''
        <RESIDUE_SELECTORS>
            <SecondaryStructure name="sheet" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="E"/>
            <SecondaryStructure name="entire_loop" overlap="0" minH="3" minE="2" include_terminal_loops="true" use_dssp="true" ss="L"/>
            <SecondaryStructure name="entire_helix" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="H"/>
            <And name="helix_cap" selectors="entire_loop">
                <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/>
            </And>
            <And name="helix_start" selectors="entire_helix">
                <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/>
            </And>
            <And name="helix" selectors="entire_helix">
                <Not selector="helix_start"/>
            </And>
            <And name="loop" selectors="entire_loop">
                <Not selector="helix_cap"/>
            </And>
        </RESIDUE_SELECTORS>
        ''')
    helix_sel = objs_sel.get_residue_selector('helix')
    loop_sel = objs_sel.get_residue_selector('loop')
    helix_cap_sel = objs_sel.get_residue_selector('helix_cap')
    
    core_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_c,helix_sel)
    bdry_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_b,helix_sel)
    surf_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_s,helix_sel)
    core_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_c,loop_sel)
    bdry_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_b,loop_sel)
    surf_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_s,loop_sel)
     
    #layer design - taskops
    core_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    core_hlx_task.aas_to_keep('AFILVW') 
    bdry_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    bdry_hlx_task.aas_to_keep('ADEHIKLNQRSTVWYM')
    surf_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    surf_hlx_task.aas_to_keep('EHKQR')
    core_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    core_loop_task.aas_to_keep('AFGILPVW') 
    bdry_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    bdry_loop_task.aas_to_keep('ADEFGHIKLNPQRSTVWY')
    surf_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    surf_loop_task.aas_to_keep('DEGHKNPQRST') 
    hlx_cap_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    hlx_cap_task.aas_to_keep('DNSTP')
    
    hlx_cap_op   = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(hlx_cap_task  , helix_cap_sel   , False)
    core_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(core_hlx_task , core_hlx_sel  , False)
    bdry_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(bdry_hlx_task , bdry_hlx_sel  , False)
    surf_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(surf_hlx_task , surf_hlx_sel  , False)
    core_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(core_loop_task, core_loop_sel, False)
    bdry_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(bdry_loop_task, bdry_loop_sel, False)
    surf_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(surf_loop_task, surf_loop_sel, False)

    #TASKFACTORY
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pack)
    tf.push_back(lock)
    tf.push_back(aroChi)
    tf.push_back(hlx_cap_op  )
    tf.push_back(core_hlx_op )
    tf.push_back(bdry_hlx_op )
    tf.push_back(surf_hlx_op )
    tf.push_back(core_loop_op)
    tf.push_back(bdry_loop_op)
    tf.push_back(surf_loop_op)
    
    
    #Design movers
    objs = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        """
        <MOVERS>
        <FastDesign name="fastdesign" repeats="1" relaxscript="KillA2019"
            cartesian="false" dualspace="false" ramp_down_constraints="false"
            bondangle="false" bondlength="false" min_type="lbfgs_armijo_nonmonotone">
        </FastDesign>
        <AddCompositionConstraintMover name="surface_polar" >
            <Comp entry="PENALTY_DEFINITION;TYPE ASP GLU HIS LYS ASN GLN ARG SER THR TYR;FRACT_DELTA_START -0.01;FRACT_DELTA_END 0.0;PENALTIES 0.1 0 ;FRACTION {};BEFORE_FUNCTION QUADRATIC;AFTER_FUNCTION CONSTANT;END_PENALTY_DEFINITION" />
        </AddCompositionConstraintMover>
        <AddCompositionConstraintMover name="ala_pen" >
                <Comp entry="PENALTY_DEFINITION;TYPE ALA;ABSOLUTE 0;PENALTIES 0 {};DELTA_START 0;DELTA_END 1;BEFORE_FUNCTION CONSTANT;AFTER_FUNCTION LINEAR;END_PENALTY_DEFINITION;" />
            </AddCompositionConstraintMover>
        </MOVERS>
        """.format(np_penalty,ala_penalty))
    surfpol = objs.get_mover('surface_polar')
    surfpol.add_residue_selector(surf_sel)
    surfpol.apply(despose)
    ala_pen = objs.get_mover('ala_pen')
    ala_pen.apply(despose)
    fd = objs.get_mover('fastdesign')
    fd.set_scorefxn(sf)
    fd.set_task_factory(tf)
    if len(diff) > 0:
        pyrosetta.rosetta.protocols.protein_interface_design.FavorNativeResidue(despose,weight)
        fd.apply(despose)
        score_list.append(sf_clean(despose))

def msd_fnr(despose,refpose,weight=1,strict_layers = 0,neighbors = 0):
    true_sel = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    allres = pyrosetta.rosetta.core.select.get_residues_from_subset(true_sel.apply(despose))
    diff = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for i in allres:
        if despose.sequence(i,i) == "C": # maintain disulfides in despose
            continue
        elif refpose.sequence(i,i) == "C": # safely replace despose residue with CYS (not CYD)
            mut = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
            mut.set_target(i)
            mut.set_res_name(pyrosetta.rosetta.core.chemical.AA(2)) # 2 is CYS
            mut.apply(despose)
        elif despose.sequence(i,i) != refpose.sequence(i,i):
            diff.append(i)
            despose.replace_residue(i,refpose.residue(i),1)
        else:
            pass
    num_des_res.append(len(diff))
    des_res.append(diff)
    #print(f"designing {len(diff)} residues")
    #what to design
    if neighbors == 1:
        designable = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(diff),6,True)
    else:
        designable = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(diff)
    packable = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector(designable,6,True)
    pack_option = pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT()
    pack = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pack_option, designable, True)
    lock_option = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    lock = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(lock_option, packable, True)
    #some task ops
    aroChi = pyrosetta.rosetta.protocols.task_operations.LimitAromaChi2Operation()
    aroChi.chi2max(110)
    aroChi.chi2min(70)
    aroChi.include_trp(True)
    #layer design - selectors -- 
    ss1 = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose1)
    ss2 = pyrosetta.rosetta.core.scoring.dssp.Dssp(pose2)
    surf_sel = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    surf_sel.set_layers(0,0,1)
    surf_sel.set_use_sc_neighbors(0)
    surf_sel.set_cutoffs(20,50)
    surf1 = pyrosetta.rosetta.core.select.get_residues_from_subset(surf_sel.apply(pose1))
    surf2 = pyrosetta.rosetta.core.select.get_residues_from_subset(surf_sel.apply(pose2))
    core_sel = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
    core_sel.set_layers(1,0,0)
    core_sel.set_use_sc_neighbors(0)
    core1 = pyrosetta.rosetta.core.select.get_residues_from_subset(core_sel.apply(pose1))
    core2 = pyrosetta.rosetta.core.select.get_residues_from_subset(core_sel.apply(pose2))
    core_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    surf_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_core = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_surf = pyrosetta.rosetta.utility.vector1_unsigned_long()
    surf_core = pyrosetta.rosetta.utility.vector1_unsigned_long()
    bdry_both = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for i in allres:
        if i in core1:
            if i in core2:
                core_both.append(i)
            elif i in surf2:
                surf_core.append(i)
            else:
                bdry_core.append(i)
        elif i in surf1:
            if i in surf2:
                surf_both.append(i)
            elif i in core2:
                surf_core.append(i)
            else:
                bdry_surf.append(i)
        else:
            if i in core2:
                bdry_core.append(i)
            elif i in surf2:
                bdry_surf.append(i)
            else:
                bdry_both.append(i)
    if len(core_both) > 0:
        sel_core_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(core_both)
    else:
        sel_core_both = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    sel_surf_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(surf_both)
    if len(bdry_core) > 0:
        sel_bdry_core = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_core)
    else:
        sel_bdry_core = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    if len(bdry_surf) > 0:
        sel_bdry_surf = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_surf)
    else:
        sel_bdry_surf = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    if len(surf_core) > 0:
        sel_surf_core = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(surf_core)
    else:
        sel_surf_core = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    sel_bdry_both = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(bdry_both)
    if strict_layers == 1:
        sel_c = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_core_both,sel_bdry_core)
        sel_b = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_bdry_both,sel_surf_core)
        sel_s = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_surf_both,sel_bdry_surf)
    else:
        sel_c = sel_core_both
        sel_s = sel_surf_both
        sel_c_or_s = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(sel_core_both,sel_surf_both)
        sel_b = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(sel_c_or_s)
        
    objs_sel = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        '''
        <RESIDUE_SELECTORS>
            <SecondaryStructure name="sheet" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="E"/>
            <SecondaryStructure name="entire_loop" overlap="0" minH="3" minE="2" include_terminal_loops="true" use_dssp="true" ss="L"/>
            <SecondaryStructure name="entire_helix" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="H"/>
            <And name="helix_cap" selectors="entire_loop">
                <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/>
            </And>
            <And name="helix_start" selectors="entire_helix">
                <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/>
            </And>
            <And name="helix" selectors="entire_helix">
                <Not selector="helix_start"/>
            </And>
            <And name="loop" selectors="entire_loop">
                <Not selector="helix_cap"/>
            </And>
        </RESIDUE_SELECTORS>
        ''')
    helix_sel = objs_sel.get_residue_selector('helix')
    loop_sel = objs_sel.get_residue_selector('loop')
    helix_cap_sel = objs_sel.get_residue_selector('helix_cap')
    
    core_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_c,helix_sel)
    bdry_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_b,helix_sel)
    surf_hlx_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_s,helix_sel)
    core_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_c,loop_sel)
    bdry_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_b,loop_sel)
    surf_loop_sel = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(sel_s,loop_sel)
     
    #layer design - taskops
    core_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    core_hlx_task.aas_to_keep('AFILVW') 
    bdry_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    bdry_hlx_task.aas_to_keep('ADEHIKLNQRSTVWYM')
    surf_hlx_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    surf_hlx_task.aas_to_keep('EHKQR')
    core_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    core_loop_task.aas_to_keep('AFGILPVW') 
    bdry_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    bdry_loop_task.aas_to_keep('ADEFGHIKLNPQRSTVWY')
    surf_loop_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    surf_loop_task.aas_to_keep('DEGHKNPQRST') 
    hlx_cap_task = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASExceptNativeRLT()
    hlx_cap_task.aas_to_keep('DNSTP')
    
    hlx_cap_op   = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(hlx_cap_task  , helix_cap_sel   , False)
    core_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(core_hlx_task , core_hlx_sel  , False)
    bdry_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(bdry_hlx_task , bdry_hlx_sel  , False)
    surf_hlx_op  = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(surf_hlx_task , surf_hlx_sel  , False)
    core_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(core_loop_task, core_loop_sel, False)
    bdry_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(bdry_loop_task, bdry_loop_sel, False)
    surf_loop_op = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(surf_loop_task, surf_loop_sel, False)

    #TASKFACTORY
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pack)
    tf.push_back(lock)
    tf.push_back(aroChi)
    tf.push_back(hlx_cap_op  )
    tf.push_back(core_hlx_op )
    tf.push_back(bdry_hlx_op )
    tf.push_back(surf_hlx_op )
    tf.push_back(core_loop_op)
    tf.push_back(bdry_loop_op)
    tf.push_back(surf_loop_op)
    
    
    #Design movers
    objs = pyrosetta.rosetta.protocols.rosetta_scripts.XmlObjects.create_from_string(
        """
        <MOVERS>
        <FastDesign name="fastdesign" repeats="1" relaxscript="KillA2019"
            cartesian="false" dualspace="false" ramp_down_constraints="false"
            bondangle="false" bondlength="false" min_type="lbfgs_armijo_nonmonotone">
        </FastDesign>
        <AddCompositionConstraintMover name="surface_polar" >
            <Comp entry="PENALTY_DEFINITION;TYPE ASP GLU HIS LYS ASN GLN ARG SER THR TYR;FRACT_DELTA_START -0.01;FRACT_DELTA_END 0.0;PENALTIES 0.1 0 ;FRACTION {};BEFORE_FUNCTION QUADRATIC;AFTER_FUNCTION CONSTANT;END_PENALTY_DEFINITION" />
        </AddCompositionConstraintMover>
        <AddCompositionConstraintMover name="ala_pen" >
                <Comp entry="PENALTY_DEFINITION;TYPE ALA;ABSOLUTE 0;PENALTIES 0 {};DELTA_START 0;DELTA_END 1;BEFORE_FUNCTION CONSTANT;AFTER_FUNCTION LINEAR;END_PENALTY_DEFINITION;" />
            </AddCompositionConstraintMover>
        </MOVERS>
        """.format(np_penalty,ala_penalty))
    surfpol = objs.get_mover('surface_polar')
    surfpol.add_residue_selector(surf_sel)
    surfpol.apply(despose)
    ala_pen = objs.get_mover('ala_pen')
    ala_pen.apply(despose)
    fd = objs.get_mover('fastdesign')
    fd.set_scorefxn(sf)
    fd.set_task_factory(tf)
    if len(diff) > 0:
        pyrosetta.rosetta.protocols.protein_interface_design.FavorNativeResidue(despose,weight)
        fd.apply(despose)
        score_list.append(sf_clean(despose))


des_y_initial(pose2,pose1,0,0)
msd_fnr(pose1,pose2,0,1,1)
msd_fnr(pose2,pose1,0,1)
msd_fnr(pose1,pose2,0.2,1)
msd_fnr(pose2,pose1,0.2,1)
msd_fnr(pose1,pose2,0.5,1)
msd_fnr(pose2,pose1,0.5,1)
msd_fnr(pose1,pose2,1,1)
msd_fnr(pose2,pose1,1,1)
msd_fnr(pose1,pose2,1.5,0)
msd_fnr(pose2,pose1,1.5,0)
msd_fnr(pose1,pose2,2,0)
msd_fnr(pose2,pose1,2,0)
np_penalty = 0.5
msd_fnr(pose1,pose2,10,0)
msd_fnr(pose2,pose1,10,0)
pose1.dump_pdb(f"X_fnr_{index}_{str(sys.argv[2])[2:]}.pdb")
pose2.dump_pdb(f"Y_fnr_{index}_{str(sys.argv[2])[2:]}.pdb")
with open("fnr_stats.list", "a+") as file:
    file.write(f"run {index}_{str(sys.argv[2])[2:]}: designed {num_des_res} residues, scores:{score_list} run4 last step 10 weight\n")

