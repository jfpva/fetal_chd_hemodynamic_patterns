# analyze_fetal_hemodynamics

Analyisis of fetal blood flow distribution patterns and blood oxygen saturations, measured using MRI, in late gestation human fetuses with a spectrum of congenital heart disease (CHD) subtypes presenting with neonatal cyanosis and age-matched fetuses with normal hearts.

`script_prep_data.m` – Matlab script to prepare data

`script_analyze_fetal_chd_flows.m` – Matlab script to run the flow analysis

`statistics.pzfx` – Prism 8 statisitcal analysis

`./results` – output of the analysis  

`./data/flow.csv` – measured flow data  

`./templates/results_template.xlsx` – Excel template for output of analysis


## 1 Disease and Anatomy Nomenclature

### 1.1 Congenital Heart Disease Sub-Types

* normal (Normal)
* hypoplastic left heart syndrome (HLHS) 
    * HLHS with restrictive atrial septum (HLHS RAS)
    * HLHS with mitral atresia and aortic stenosis (HLHS MA AS)
    * HLHS with mitral stenosis and aortic stenosis (HLHS MS AS)
    * HLHS with mitral stenosis and aortic atresia (HLHS MS AA)
    * HLHS with mitral atresia and aortic atresia (HLHS MA AA)
    * HLHS with double outlet right ventricle (HLHS DORV)
* transposition of the great arteries (TGA)
    * TGA with intact ventricular septum (TGA IVS)
    * TGA with ventricular septal defect (TGA VSD)
    * TGA with VSD and pulmonary stenosis (TGA VSD PS)
    * TGA with intact ventricular septum with coarctation of the aorta (TGA COA)
* tetralogy of Fallot (TOF)
    * with pulmonary stenosis (TOF ONLY)
    * with pulmonary atresia (TOF PA)
* Ebstein's anomaly (EA)
    * EA with functional pulmonary atresia with no circular shunt (Ebstein's no Circular Shunt)
    * EA with a circular shunt (Ebstein's Circular Shunt)
* tricuspid atresia (TA)
    * TA with ventriculo-arterial concordance (TA VA Concordance)
    * TA with ventriculo-arterial discordance (TA VA Discordance)

### 1.2 Fetal Vessels

* __AAo__ – ascending aorta
* __DAo__ – desending aorta
* __DA__ – ductus arteriosus
* __MPA__ – main pulmonary artery
* __PBF__ – pulmonary blood flow (sum of left and right pulmonary arteries)
* __SVC__ – superior vena cava
* __IVC__ – inferior vena cava
* __UV__ – umbilical vein
* __CA__ – coronary artery
* __CS__ – coronary sinus
* __FO__ – foramen ovale
* __ICS__ – intracardiac shunt
* __CVO__ – combined ventricular output


## 2 Derived Flows

### 2.1 Missing Measurements

Some measurements may be missing in particularly challenging cases. 
In these circumstances, it may be possible to derive missing measurements from measured flows using conservation of flow at the arch junciton and/or pulmonary branch junction in order to calculate CVO.

__Arch Junction:__ At the junction of the aorta and the DA, `Q_AAo = Q_DAo - Q_DA + Q_SVC`.

__Pulmonary Branch Junction:__ At the junction of the MPA and branch pulmonary arteries, `Q_MPA = Q_DA + Q_PBF`.

These relationships were used to derive missing measurements, first, directly from measured flows and second, using other derived flow.

#### Outliers

Any derived flows which produced retrograde flow in the MPA, AAo, DAo, SVC, PBF or UV were considered outliers and discarded, except for MPA flow in cases with Ebstein's anomaly with circular shunt and AAo flow in cases with aortic atresia. 
In cases with Ebstein's anomaly or pulmonary atresia, derived antegrade DA flows were also considered outliers.

Furthermore, in sub-groups with sufficient number of measurements to estimate population distributions, any derived flows that were more than three scaled median absolute deviations from the median where considered outliers. Details of this approach can be found [here](https://www.mathworks.com/help/matlab/ref/isoutlier.html).

### 2.2 Combined Ventricular Output

Combined ventricular output (CVO) can be used to interpret the distribution of blood flow as a fraction of the total blood flow. 

#### Normal Circulation

In normal circulation, `Q_CVO = Q_MPA + Q_AAo + Q_CA` and can be calculated as `Q_CVO = ( Q_MPA + Q_AAo ) / 0.97` where coronary blood flow is estimated as 3% of CVO [Rudulolph2001]. 
Hence, coronary artery (CA) blood flow can be calculated as `Q_{CA} = 0.03 Q_CVO`, which is equal to coronary sinus (CS) flow.

This calculation holds for all circulation types with antegrade MPA and AAo flow.

#### Pulmonary Atresia

In circulations with pulmonary atresia, `Q_CVO = Q_AAo / 0.97` and `Q_MPA = 0`.

#### Aortic Atresia

In circulations with aortic atresia `Q_CVO = Q_MPA` and the coronary areteries are supplied by retrograde flow through the ascending aorta, i.e., `Q_AAo = - Q_{CA}.`

### 2.3 Inferior Vena Cava Flow

`Q_IVC = Q_DAo`  

### 2.4 Foramen Ovale Blood Flow

Calculation of FO flow is dependent on the type of circulation.

__Normal__  
Normal, HLHS MS AS  
`Q_FO = Q_AAo + Q_CA - Q_PBF`  

__Transposition__  
TGA IVS, TGA COA  
`Q_FO = Q_SVC + Q_IVC - Q_AAo`
            
__Aortic Atresia__  
HLHS RAS, HLHS MS AA, HLHS MA AA  
`Q_FO = - Q_PBF`  

__Tricuspid Atresia, Functional Pulmonary Atresia__ 
TA VA, Ebstein's Anomaly with no Circular Shunt  
`Q_FO = Q_SVC + Q_IVC + Q_CS`  

__Ebstein's Anomaly with Circular Shunt__
`Q_FO = -Q_MPA + Q_SVC + Q_IVC + Q_CS`

__Ventricular Septal Defect, Overriding Aorta__  
HLHS MA AS, HLHS DORV, TGA VSD, TGA VSD PS, TOF  
`Q_FO` unknown


## 3 Flow Distribution Models


Models of the distribution of blood flow, as a percentage of CVO, were extrapolated from measured flows based on a constrained nonlinear optimization satisfying conservation of flow throughout the fetal circulatory system using a weighted root mean squared difference objective function to limit the change from measured values. 

`F_objective(x) = ( sum( w * ( x - x0 )^2 ) / sum( w ) )^0.5`, where `x0` are the median measured flows as %CVO for a group or sub-group, `x` are distribution model flows that satisfy the constraints and `w` are weights. 

Equality constraints were based on conservation of blood flow at the following locations:

1. blood flow in and out of the entire heart
2. pulmonary artery branch junction
3. arch junction
4. right ventriculo-arterial connection
5. left ventriculo-arterial connection
6. lower body and placenta circulation
7. coronary circulation
8. ventricular output
9. atrial input
10. combined ventricular output
11. combined atrial input

Specifics of equality and inequality constraints, as well as bounds and weights were as follows.

Intracardiac shunt (ICS) was used to quantify overall right to left shunt, which is equivalent to FO flow in circulations without ventricular septal defect or overriding aorta.

### circulation type: normal

```
Equality Constraints:
 - Q_mpa - Q_aao + Q_svc + Q_pbf + Q_ivc - Q_ca  + Q_cs        =            0     
 + Q_mpa - Q_da  - Q_pbf                                       =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 - Q_mpa + Q_svc - Q_ics + Q_ivc + Q_cs                        =            0     
 - Q_aao + Q_pbf + Q_ics - Q_ca                                =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_mpa + Q_aao - 32.3333 Q_ca                                =            0     
 + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                        =            0     
 + Q_mpa + Q_aao + Q_ca                                        =          100     
 + Q_svc + Q_pbf + Q_ivc + Q_cs                                =          100     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower       0     0     0  -100     0     0     0  -100     0     0     0
    upper     100   100   100   100   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.60  0.60  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00
```


### circulation type: transposition

```
Equality Constraints:
 - Q_mpa - Q_aao + Q_svc + Q_pbf + Q_ivc - Q_ca  + Q_cs        =            0     
 + Q_mpa - Q_da  - Q_pbf                                       =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 - Q_aao + Q_svc - Q_ics + Q_ivc - Q_ca  + Q_cs                =            0     
 - Q_mpa + Q_pbf + Q_ics                                       =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_mpa + Q_aao - 32.3333 Q_ca                                =            0     
 + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                        =            0     
 + Q_mpa + Q_aao + Q_ca                                        =          100     
 + Q_svc + Q_pbf + Q_ivc + Q_cs                                =          100     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower       0     0     0  -100     0     0     0  -100     0     0     0
    upper     100   100   100   100   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.60  0.60  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00
```

### circulation type: hlhs_aortic_atresia

```
Equality Constraints:
 - Q_mpa + Q_svc + Q_pbf + Q_ivc + Q_cs                        =            0     
 + Q_mpa - Q_da  - Q_pbf                                       =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 - Q_mpa + Q_svc - Q_ics + Q_ivc + Q_cs                        =            0     
 + Q_pbf + Q_ics                                               =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                        =            0     
 + Q_mpa                                                       =          100     
 + Q_svc + Q_pbf + Q_ivc + Q_cs                                =          100     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower       0  -100     0  -100     0     0     0  -100     0     0     0
    upper     100     0   100   100   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.60  0.60  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00         
```

### circulation type: tof_pa

```
Equality Constraints:
 - Q_aao + Q_svc + Q_pbf + Q_ivc - Q_ca  + Q_cs                =            0     
 - Q_da  - Q_pbf                                               =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 + Q_svc - Q_ics + Q_ivc + Q_cs                                =            0     
 - Q_aao + Q_pbf + Q_ics - Q_ca                                =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_aao - 32.3333 Q_ca                                        =            0     
 + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                        =            0     
 + Q_aao + Q_ca                                                =          100     
 + Q_svc + Q_pbf + Q_ivc + Q_cs                                =          100     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower       0     0     0  -100     0     0     0  -100     0     0     0
    upper       0   100   100     0   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.00  0.60  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00
```

### circulation type: double_outlet_right_ventricle

```
Equality Constraints:
 - Q_mpa - Q_aao + Q_svc + Q_pbf + Q_ivc - Q_ca  + Q_cs        =            0     
 + Q_mpa - Q_da  - Q_pbf                                       =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 - Q_mpa - Q_aao + Q_svc - Q_ics + Q_ivc - Q_ca  + Q_cs        =            0     
 + Q_pbf + Q_ics                                               =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_mpa + Q_aao - 32.3333 Q_ca                                =            0     
 + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                        =            0     
 + Q_mpa + Q_aao + Q_ca                                        =          100     
 + Q_svc + Q_pbf + Q_ivc + Q_cs                                =          100     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower       0     0     0  -100     0     0     0  -100     0     0     0
    upper     100   100   100   100   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.60  0.30  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00
```

### circulation type: ebsteins_circular_shunt

```
Equality Constraints:
 - Q_mpa - Q_aao + Q_svc + Q_pbf + Q_ivc - Q_ca  + Q_cs        =            0     
 + Q_mpa - Q_da  - Q_pbf                                       =            0     
 + Q_aao - Q_svc + Q_da  - Q_dao                               =            0     
 - Q_mpa + Q_svc - Q_ics + Q_ivc + Q_cs                        =            0     
 - Q_aao + Q_pbf + Q_ics - Q_ca                                =            0     
 + Q_dao - Q_ivc                                               =            0     
 + Q_ca  - Q_cs                                                =            0     
 + Q_aao - 32.3333 Q_ca                                        =            0     
 - Q_mpa + Q_svc + Q_pbf + Q_ivc - 32.3333 Q_cs                =            0     

Inequality Constraints:
 - Q_dao + Q_uv                                               <=            0

Bounds:       mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
    lower    -100     0     0  -100     0     0     0  -100     0     0     0
    upper       0   100   100     0   100   100   100   100   100   100   100

Weights:      mpa   aao   svc    da   dao   pbf    uv   ics   ivc    ca    cs
             0.30  0.60  0.80  0.50  1.00  0.28  1.00  0.00  0.00  0.00  0.00
```
