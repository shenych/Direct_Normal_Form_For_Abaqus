##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################




from abaqus import *
from abaqusConstants import *
import time




session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=170.88020324707, 
    height=134.375)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)



Inputpath=Current_path+"inp_files/"
Outputpath=Current_path+"rpt_files/"
Odbpath=Current_path



text='("RF", NODAL, ('

for i in range(1,dof_RF+1):
    text=text+'(COMPONENT, "RF'+str(i)+'"),'

text=text+')),'

if dof_RM>0:
    text=text+'("RM", NODAL, ('
    for j in range(1,dof_RM+1):
        text=text+'(COMPONENT, "RM'+str(j)+'"),'
    text=text+')),'



tensortype='absd'
for loop in [0,1]:

    for mode_i in range(Nmodes):
        for mode_j in range(mode_i,Nmodes):    
            for ne_or_pi in [1,2]:
                Inputname=tensortype[loop]+'_ten_'+str(mode_i+1)+str(mode_j+1)+'_'+str(ne_or_pi) 
                mdb.ModelFromInputFile(name=Inputname, 
                               inputFileName=Inputpath+Inputname+'.inp')

                session.viewports['Viewport: 1'].assemblyDisplay.setValues(
                        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
                a = mdb.models[Inputname].rootAssembly
                session.viewports['Viewport: 1'].setValues(displayedObject=a)
                mdb.Job(name=Inputname, model=Inputname, description='', 
                        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                        explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
                        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
                        multiprocessingMode=MPI, numCpus=4, numDomains=4,numGPUs=0)
                mdb.jobs[Inputname].submit(consistencyChecking=OFF)



  
 

    for mode in range(Nmodes):
        for mode_i in range(Nmodes):
            for mode_j in range(mode_i,Nmodes):  
                for ne_or_pi in [1,2]:
                    Inputname='Phi_'+str(mode+1)+'_'+tensortype[loop]+'_ten_'+str(mode_i+1)+str(mode_j+1)+'_'+str(ne_or_pi)   
                    mdb.ModelFromInputFile(name=Inputname, 
                               inputFileName=Inputpath+Inputname+'.inp')

                    session.viewports['Viewport: 1'].assemblyDisplay.setValues(
                        optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
                    a = mdb.models[Inputname].rootAssembly
                    session.viewports['Viewport: 1'].setValues(displayedObject=a)
                    mdb.Job(name=Inputname, model=Inputname, description='', 
                        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
                        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
                        explicitPrecision=SINGLE, nodalOutputPrecision=FULL, echoPrint=OFF, 
                        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
                        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
                        multiprocessingMode=MPI, numCpus=4, numDomains=4,numGPUs=0)
                    mdb.jobs[Inputname].submit(consistencyChecking=OFF)


time.sleep(15)


for loop in [0,1]:
    for mode_i in range(Nmodes):
        for mode_j in range(mode_i,Nmodes):    
            for ne_or_pi in [1,2]:
                Inputname=tensortype[loop]+'_ten_'+str(mode_i+1)+str(mode_j+1)+'_'+str(ne_or_pi) 


                o3 = session.openOdb(name=Odbpath+Inputname+'.odb')

                session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                session.viewports['Viewport: 1'].makeCurrent()
                odb = session.odbs[Odbpath+Inputname+'.odb']
                nf = NumberFormat(numDigits=9, precision=0, format=ENGINEERING)
                session.fieldReportOptions.setValues(printTotal=OFF, printMinMax=OFF, numberFormat=nf)
                session.writeFieldReport(
                    fileName=Outputpath+Inputname+'.rpt', 
                    append=ON, sortItem='Node Label', odb=odb, step=0, frame=10, 
                    outputPosition=NODAL, variable=(eval(text)), stepFrame=SPECIFY)


  
 

    for mode in range(Nmodes):
        for mode_i in range(Nmodes):
            for mode_j in range(mode_i,Nmodes):  
                for ne_or_pi in [1,2]:
                    Inputname='Phi_'+str(mode+1)+'_'+tensortype[loop]+'_ten_'+str(mode_i+1)+str(mode_j+1)+'_'+str(ne_or_pi)   

                    o3 = session.openOdb(name=Odbpath+Inputname+'.odb')

                    session.viewports['Viewport: 1'].setValues(displayedObject=o3)
                    session.viewports['Viewport: 1'].makeCurrent()
                    odb = session.odbs[Odbpath+Inputname+'.odb']
                    nf = NumberFormat(numDigits=9, precision=0, format=ENGINEERING)
                    session.fieldReportOptions.setValues(printTotal=OFF, printMinMax=OFF, numberFormat=nf)
                    session.writeFieldReport(
                        fileName=Outputpath+Inputname+'.rpt', 
                        append=ON, sortItem='Node Label', odb=odb, step=0, frame=10, 
                        outputPosition=NODAL, variable=(eval(text)), stepFrame=SPECIFY)

 