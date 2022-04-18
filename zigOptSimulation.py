from abaqus import *
from abaqusConstants import *
import __main__

def zigOptSimulation(profile, DPC):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    import numpy as np

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=30.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(viewStyle=AXISYM)
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -15.0), point2=(0.0, 15.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(0.0, 10.0), point2=(0.0, 0.0))
    s1.VerticalConstraint(entity=g[3], addUndoState=False)
    for a, i in enumerate(profile):
        if i == 0:
            s1.Line(point1=(0.0,0.0), point2=(0.0,profile[i]))
        elif i == len(profile)-1:
            s1.Line(point1=(10.0, profile[i]), point2=(0.0,10.0))
        else:
            s1.Line(point1=(0.5*(i-1),profile[i-1]),point2=(0.5*i,profile[i]))
    s1.HorizontalConstraint(entity=g[24], addUndoState=False)
    p = mdb.models['Model-1'].Part(name='zig', dimensionality=AXISYMMETRIC, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['zig']
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    #create part: zig

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=100.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(viewStyle=AXISYM)
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -50.0), point2=(0.0, 50.0))
    s1.FixedConstraint(entity=g[2])
    s1.rectangle(point1=(0.0, 0.0), point2=(29.9, 40.9))
    p = mdb.models['Model-1'].Part(name='specimen', dimensionality=AXISYMMETRIC, 
        type=DEFORMABLE_BODY)
    p = mdb.models['Model-1'].parts['specimen']
    p.BaseShell(sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    #create part: specimen 

    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=100.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.sketchOptions.setValues(viewStyle=AXISYM)
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -50.0), point2=(0.0, 50.0))
    s1.FixedConstraint(entity=g[2])
    s1.Line(point1=(0.0, 0.0), point2=(29.9, 0.0))
    s1.HorizontalConstraint(entity=g[3], addUndoState=False)
    s1.Line(point1=(29.9, 0.0), point2=(29.9, 50.0))
    s1.VerticalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    p = mdb.models['Model-1'].Part(name='holder', dimensionality=AXISYMMETRIC, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts['holder']
    p.BaseWire(sketch=s1)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']
    #create part: holder

    mdb.models['Model-1'].Material(name='granularMaterials')
    mdb.models['Model-1'].materials['granularMaterials'].Density(table=((1.3e-09, 
        ), ))
    mdb.models['Model-1'].materials['granularMaterials'].Elastic(table=((DPC[8], 
        DPC[9]), ))
    mdb.models['Model-1'].materials['granularMaterials'].CapPlasticity(table=((
        DPC[7], DPC[6], DPC[4], 0.0, DPC[5], 1.0), ))
    def pv_strain(pb, l=DPC[0], k=DPC[1], A=DPC[2], pb0=DPC[3]):
        return (l-k)*np.log((1+A-l*np.log(1000*pb0))/(1+A-l*np.log(1000*pb)))/l
    table = ((pb,pv_strain(pb)) for pb in np.arange(DPC[3],DPC[3]+0.005*1000, 0.005))
    #create material: granularMaterials

    mdb.models['Model-1'].materials['granularMaterials'].capPlasticity.CapHardening(
        table=table)
    mdb.models['Model-1'].Material(name='steel')
    mdb.models['Model-1'].materials['steel'].Density(table=((8.05e-09, ), ))
    mdb.models['Model-1'].materials['steel'].Elastic(table=((200000.0, 0.27), ))
    #create material: steel

    mdb.models['Model-1'].HomogeneousSolidSection(name='specimenSection', 
        material='granularMaterials', thickness=None)
    #create section: sepcimenSection

    mdb.models['Model-1'].HomogeneousSolidSection(name='zig_holderSection', 
        material='steel', thickness=None)
    #create section: zig_holderSection

    p = mdb.models['Model-1'].parts['specimen']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(faces=faces, name='Set-1')
    p.SectionAssignment(region=region, sectionName='specimenSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    #assign section: specimenSection to specimen

    p = mdb.models['Model-1'].parts['zig']
    f = p.faces
    faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
    region = p.Set(faces=faces, name='Set-1')
    p.SectionAssignment(region=region, sectionName='zig_holderSection', offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    #assign section: zig_holderSection to zig

    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByThreePoints(coordSysType=CYLINDRICAL, origin=(0.0, 0.0, 0.0), 
        point1=(1.0, 0.0, 0.0), point2=(0.0, 0.0, -1.0))
    p = mdb.models['Model-1'].parts['holder']
    a.Instance(name='holder-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['specimen']
    a.Instance(name='specimen-1', part=p, dependent=ON)
    p = mdb.models['Model-1'].parts['zig']
    a.Instance(name='zig-1', part=p, dependent=ON)
    #initialize assembly

    a = mdb.models['Model-1'].rootAssembly
    a.translate(instanceList=('zig-1', ), vector=(0.0, 41, 0.0))
    #translate zig
    
    mdb.models['Model-1'].ExplicitDynamicsStep(name='pressing', previous='Initial', 
        timePeriod=12.0)
    #create step: pressing
    
    p1 = mdb.models['Model-1'].parts['zig']
    session.viewports['Viewport: 1'].setValues(displayedObject=p1)
    p = mdb.models['Model-1'].parts['zig']
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=v1[2])
    r = p.referencePoints
    refPoints=(r[3], )
    p.Set(referencePoints=refPoints, name='RP')
    #create reference point: RP in zig

    regionDef=mdb.models['Model-1'].rootAssembly.allInstances['zig-1'].sets['RP']
    mdb.models['Model-1'].historyOutputRequests['H-Output-1'].setValues(variables=(
        'U2', 'RF2'), timeInterval=0.1, region=regionDef, 
        sectionPoints=DEFAULT, rebar=EXCLUDE)
    #create history output: H-Output-1
    
    mdb.models['Model-1'].ContactProperty('fric02')
    mdb.models['Model-1'].interactionProperties['fric02'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.2, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    #create interaction property: fric02

    p = mdb.models['Model-1'].parts['zig']
    s = p.edges
    side1Edges = s.getSequenceFromMask(mask=('[#3ffffc ]', ), )
    p.Surface(side1Edges=side1Edges, name='shape')
    #create surface: shape on zig
    
    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['zig-1'].surfaces['shape']
    s1 = a.instances['specimen-1'].edges
    side1Edges1 = s1.getSequenceFromMask(mask=('[#4 ]', ), )
    region2=a.Surface(side1Edges=side1Edges1, name='s_Surf-1')
    mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='zig_specimen', 
        createStepName='pressing', master = region1, slave = region2, 
        mechanicalConstraint=PENALTY, sliding=FINITE, 
        interactionProperty='fric02', initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
    #create interaction: zig_specimen

    p = mdb.models['Model-1'].parts['specimen']
    s = p.edges
    side1Edges = s.getSequenceFromMask(mask=('[#3 ]', ), )
    p.Surface(side1Edges=side1Edges, name='holderContact')
    #create surface: holderContact on specimen
    
    p = mdb.models['Model-1'].parts['holder']
    s = p.edges
    side1Edges = s.getSequenceFromMask(mask=('[#3 ]', ), )
    p.Surface(side1Edges=side1Edges, name='whole')
    #create surface: whole on holder

    a = mdb.models['Model-1'].rootAssembly
    region1=a.instances['holder-1'].surfaces['whole']
    region2=a.instances['specimen-1'].surfaces['holderContact']
    mdb.models['Model-1'].SurfaceToSurfaceContactExp(name ='specimen_holder', 
        createStepName='pressing', master = region1, slave = region2, 
        mechanicalConstraint=PENALTY, sliding=FINITE, 
        interactionProperty='fric02', initialClearance=OMIT, datumAxis=None, 
        clearanceRegion=None)
    #create interaction: specimen_holder

    a = mdb.models['Model-1'].rootAssembly
    f1 = a.instances['zig-1'].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
    region2=a.Set(faces=faces1, name='b_Set-2')
    r1 = a.instances['zig-1'].referencePoints
    refPoints1=(r1[3], )
    region1=regionToolset.Region(referencePoints=refPoints1)
    mdb.models['Model-1'].RigidBody(name='Constraint-1', refPointRegion=region1, 
        bodyRegion=region2)
    #create rigid body constraint: Constraint-1
    
    a = mdb.models['Model-1'].rootAssembly
    e1 = a.instances['specimen-1'].edges
    edges1 = e1.getSequenceFromMask(mask=('[#8 ]', ), )
    region = a.Set(edges=edges1, name='Set-4')
    mdb.models['Model-1'].XsymmBC(name='specimenXsymm', createStepName='Initial', 
        region=region, localCsys=None)
    #create boundary condition: specimenXsymm
    
    p = mdb.models['Model-1'].parts['holder']
    v2, e1, d2, n1 = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=v2[0])
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['holder-1'].referencePoints
    refPoints1=(r1[3], )
    region = a.Set(referencePoints=refPoints1, name='Set-5')
    mdb.models['Model-1'].EncastreBC(name='holder_encastre', 
        createStepName='Initial', region=region, localCsys=None)
    #create boundary condition: holder_encastre

    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['zig-1'].referencePoints
    refPoints1=(r1[3], )
    region = a.Set(referencePoints=refPoints1, name='Set-6')
    mdb.models['Model-1'].XsymmBC(name='zigXsymm', createStepName='Initial', 
        region=region, localCsys=None)
    #create boundary condition: zigXsymm
    
    a = mdb.models['Model-1'].rootAssembly
    r1 = a.instances['zig-1'].referencePoints
    refPoints1=(r1[3], )
    region = a.Set(referencePoints=refPoints1, name='Set-7')
    mdb.models['Model-1'].VelocityBC(name='zigVelocity', createStepName='pressing', 
        region=region, v1=UNSET, v2=-0.833, vr3=UNSET, amplitude=UNSET, 
        localCsys=None, distributionType=UNIFORM, fieldName='')
    #create boundary condition: zigVelocity
    
    p = mdb.models['Model-1'].parts['holder']
    p.seedPart(size=5.0, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()
    #create mesh: holder

    p = mdb.models['Model-1'].parts['specimen']
    p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=QUAD, technique=STRUCTURED)
    elemType1 = mesh.ElemType(elemCode=CAX4R, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, hourglassControl=DEFAULT, 
        distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=CAX3, elemLibrary=STANDARD)
    pickedRegions =(pickedRegions, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()
    #create mesh: specimen
    
    p = mdb.models['Model-1'].parts['zig']
    p.seedPart(size=0.5, deviationFactor=0.1, minSizeFactor=0.1)
    f = p.faces
    pickedRegions = f.getSequenceFromMask(mask=('[#1 ]', ), )
    p.setMeshControls(regions=pickedRegions, elemShape=TRI)
    elemType1 = mesh.ElemType(elemCode=CAX4R, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=CAX3, elemLibrary=STANDARD, 
        secondOrderAccuracy=OFF, distortionControl=DEFAULT)
    pickedRegions =(pickedRegions, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
    p.generateMesh()
    #create mesh: zig

    mdb.Job(name='pressing', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK, 
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=2, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=2)
    #create job: pressing

    mdb.jobs['pressing'].submit(consistencyChecking=OFF)
    #submit job

