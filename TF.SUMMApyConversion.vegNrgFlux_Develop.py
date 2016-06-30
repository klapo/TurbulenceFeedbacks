# SUMMA - Structure for Unifying Multiple Modeling Alternatives
# Copyright (C) 2014-2015 NCAR/RAL
#
# This file is part of SUMMA
#
# For more information see: http://www.ral.ucar.edu/projects/summa
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# *******************************************************************************************************
# python conversion of SUMMA's subroutine vegNrgFlux. Vegetation/canopy is neglected. 
# *******************************************************************************************************

## Define variables

# input: model control
#firstSubStep,                            # intent(in): flag to indicate if we are processing the first sub-step
#firstFluxCall,                           # intent(in): flag to indicate if we are processing the first flux call
#
## input: model state variables
# upperBoundTemp,                          # intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
# groundTempTrial,                         # intent(in): trial value of ground temperature (K)
#
## input: model derivatives
#
## input/output: data structures
# forc_data,                               # intent(in):    model forcing data
# mpar_data,                               # intent(in):    model parameters
# mvar_data,                               # intent(inout): model variables for a local HRU
# bvar_data,                               # intent(in):    model variables for the local basin
# model_decisions,                         # intent(in):    model decisions
#
## output: liquid water fluxes associated with evaporation/transpiration (needed for coupling)
# returnGroundEvaporation,                 # intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
#
## output: fluxes
# groundNetFlux,                           # intent(out): net energy flux for the ground surface (W m-2)
#
## output: energy flux derivatives
# dGroundNetFlux_dGroundTemp,              # intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
#
# # input: model decisions
# ix_bcUpprTdyn                     # intent(in): [i4b] choice of upper boundary condition for thermodynamics
# ix_fDerivMeth                     # intent(in): [i4b] choice of method to compute derivatives
# ix_astability                     # intent(in): [i4b] choice of stability function
# ix_groundwatr                     # intent(in): [i4b] choice of groundwater parameterization
#
# # input: physical attributes
# soilTypeIndex                     # intent(in): [i4b] soil type index
# # input: aerodynamic resistance parameters
# z0Snow                          # intent(in): [dp] roughness length of snow (m)
# z0Soil                          # intent(in): [dp] roughness length of soil (m)
# critRichNumber                  # intent(in): [dp] critical value for the bulk Richardson number where turbulence ceases (-)
# Louis79_bparam                  # intent(in): [dp] parameter in Louis (1979) stability function
# Louis79_cStar                   # intent(in): [dp] parameter in Louis (1979) stability function
# Mahrt87_eScale                  # intent(in): [dp] exponential scaling factor in the Mahrt (1987) stability function
#
# # input: soil stress parameters
# theta_sat                       # intent(in): [dp] soil porosity (-)
# theta_res                       # intent(in): [dp] residual volumetric liquid water content (-)
# plantWiltPsi                    # intent(in): [dp] matric head at wilting point (m)
# soilStressParam                 # intent(in): [dp] parameter in the exponential soil stress function (-)
# critSoilWilting                 # intent(in): [dp] critical vol. liq. water content when plants are wilting (-)
# critSoilTranspire               # intent(in): [dp] critical vol. liq. water content when transpiration is limited (-)
# critAquiferTranspire            # intent(in): [dp] critical aquifer storage value when transpiration is limited (m)
# minStomatalResistance           # intent(in): [dp] mimimum stomatal resistance (s m-1)
#
# # input: forcing at the upper boundary
# mHeight                         # intent(in): [dp] measurement height (m)
# airtemp                         # intent(in): [dp] air temperature at some height above the surface (K)
# windspd                         # intent(in): [dp] wind speed at some height above the surface (m s-1)
# airpres                         # intent(in): [dp] air pressure at some height above the surface (Pa)
# LWRadAtm                        # intent(in): [dp] downwelling longwave radiation at the upper boundary (W m-2)
# scalarVPair                     # intent(in): [dp] vapor pressure at some height above the surface (Pa)
# scalarTwetbulb                  # intent(in): [dp] wetbulb temperature (K)
# scalarRainfall                  # intent(in): [dp] computed rainfall rate (kg m-2 s-1)
# scalarSnowfall                  # intent(in): [dp] computed snowfall rate (kg m-2 s-1)
#
# # input: water storage
# # NOTE: soil stress only computed at the start of the substep (firstFluxCall=.true.)
# scalarSWE                       # intent(in): [dp]    snow water equivalent on the ground (kg m-2)
# scalarSnowDepth                 # intent(in): [dp]    snow depth on the ground surface (m)
# mLayerVolFracLiq                # intent(in): [dp(:)] volumetric fraction of liquid water in each soil layer (-)
# mLayerMatricHead                # intent(in): [dp(:)] matric head in each layer (m)
# localAquiferStorage             # intent(in): [dp]    aquifer storage for the local column (m)
# basinAquiferStorage             # intent(in): [dp]    aquifer storage for the single basin (m)
#
# # input: shortwave radiation fluxes
# scalarGroundAbsorbedSolar       # intent(in): [dp] solar radiation absorbed by ground (W m-2)
#
# # output: fraction of wetted canopy area and fraction of snow on the ground
# scalarGroundSnowFraction        # intent(out): [dp] fraction of ground covered with snow (-)
#
# # output: longwave radiation fluxes
# scalarLWRadGround               # intent(out): [dp] longwave radiation emitted at the ground surface (W m-2)
# scalarLWRadUbound2Ground        # intent(out): [dp] downward atmospheric longwave radiation absorbed by the ground (W m-2)
# scalarLWRadUbound2Ubound        # intent(out): [dp] atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
# scalarLWRadGround2Ubound        # intent(out): [dp] longwave radiation emitted from ground lost thru upper boundary (W m-2)
# scalarLWNetGround               # intent(out): [dp] net longwave radiation at the ground surface (W m-2)
# scalarLWNetUbound               # intent(out): [dp] net longwave radiation at the upper boundary (W m-2)
#
# # output: aerodynamic resistance
# scalarRiBulkGround              # intent(out): [dp] bulk Richardson number for the ground surface (-)
# scalarFrictionVelocity          # intent(out): [dp] friction velocity (m s-1)
# scalarGroundResistance          # intent(out): [dp] below canopy aerodynamic resistance (s m-1)
#
# # input/output: soil resistance  call
# mLayerRootDensity               # intent(in):    [dp] root density in each layer (-)
# scalarAquiferRootFrac           # intent(in):    [dp] fraction of roots below the lowest soil layer (-)
# scalarTranspireLim              # intent(inout): [dp] weighted average of the transpiration limiting factor (-)
# mLayerTranspireLim              # intent(inout): [dp] transpiration limiting factor in each layer (-)
# scalarTranspireLimAqfr          # intent(inout): [dp] transpiration limiting factor for the aquifer (-)
# scalarSoilRelHumidity           # intent(inout): [dp] relative humidity in the soil pores [0-1]
# scalarSoilResistance            # intent(inout): [dp] resistance from the soil (s m-1)
#
# # output: turbulent heat fluxes
# scalarLatHeatSubVapGround       # intent(inout): [dp] latent heat of sublimation/vaporization for the ground surface (J kg-1)
# scalarSatVP_groundTemp          # intent(out):   [dp] saturation vapor pressure at the temperature of the ground surface (Pa)
# scalarSenHeatTotal              # intent(out):   [dp] sensible heat from the canopy air space to the atmosphere (W m-2)
# scalarSenHeatGround             # intent(out):   [dp] sensible heat flux from ground surface below vegetation (W m-2)
# scalarLatHeatTotal              # intent(out):   [dp] latent heat from the canopy air space to the atmosphere (W m-2)
# scalarLatHeatGround             # intent(out):   [dp] latent heat flux from ground surface below vegetation (W m-2)
#
# # output: advective heat fluxes
# scalarCanopyAdvectiveHeatFlux   # intent(out): [dp] heat advected to the canopy surface with rain + snow (W m-2)
# scalarGroundAdvectiveHeatFlux   # intent(out): [dp] heat advected to the ground surface with throughfall (W m-2)
#
# # output: mass fluxes
# scalarCanopySublimation         # intent(out): [dp] canopy sublimation/frost (kg m-2 s-1)
# scalarSnowSublimation           # intent(out): [dp] snow sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)
#
# # input/output: canopy air space
# scalarVP_CanopyAir              # intent(inout): [dp] vapor pressure of the canopy air space (Pa)
# scalarCanopyStabilityCorrection # intent(inout): [dp] stability correction for the canopy (-)
# scalarGroundStabilityCorrection # intent(inout): [dp] stability correction for the ground surface (-)
#
# # output: liquid water fluxes
# scalarCanopyTranspiration       # intent(out): [dp] canopy transpiration (kg m-2 s-1)
# scalarCanopyEvaporation         # intent(out): [dp] canopy evaporation/condensation (kg m-2 s-1)
# scalarGroundEvaporation         # intent(out): [dp] ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)


# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# ***** PRELIMINARIES  **********************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************

# flux boundary condition
case(energyFlux)

# identify the appropriate groundwater variable
scalarGroundStabilityCorrection_old = scalarGroundStabilityCorrection       # stability correction for the ground surface (-)

# set latent heat of sublimation/vaporization for canopy and ground surface (Pa/K)
if firstFluxCall:

	# case when there is snow on the ground (EXCLUDE "snow without a layer" -- in this case, evaporate from the soil)
	if nSnow > 0:
		if groundTempTrial > Tfreeze:
			error('do not expect ground temperature > 0 when snow is on the ground')
		scalarLatHeatSubVapGround = LH_sub  # sublimation from snow
		scalarGroundSnowFraction  = 1.
	# case when the ground is snow-free
	else:
		scalarLatHeatSubVapGround = LH_vap  # evaporation of water in the soil pores: this occurs even if frozen because of super-cooled water
		scalarGroundSnowFraction  = 0.

# compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
z0Ground = z0soil*(1. - scalarGroundSnowFraction) + z0Snow*scalarGroundSnowFraction     # roughness length (m)

# compute emissivity of the ground surface (-)
groundEmissivity = scalarGroundSnowFraction*snowEmissivity + (1. - scalarGroundSnowFraction)*soilEmissivity  # emissivity of the ground surface (-)

# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# ***** AERODYNAMIC RESISTANCE *****************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************

# NOTE: compute for all iterations

# compute aerodynamic resistances
# Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
#       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
#       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
aeroOut = aeroResist(\
	# input: model control
	(ix_fDerivMeth == analytical), \     # logical flag if would like to compute analytical derivaties
	ix_astability,                 \     # choice of stability function
	# input: above-canopy forcing data
	mHeight,                       \     # measurement height (m)
	airtemp,                       \     # air temperature at some height above the surface (K)
	windspd,                       \     # wind speed at some height above the surface (m s-1)
	# input: canopy and ground temperature
	groundTempTrial,               \     # temperature of the ground surface (K)
	# input: diagnostic variables
	scalarSnowDepth,               \     # snow depth (m)
	# input: parameters
	z0Ground,                      \     # roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
	critRichNumber,                \     # critical value for the bulk Richardson number where turbulence ceases (-)
	Louis79_bparam,                \     # parameter in Louis (1979) stability function
	Mahrt87_eScale,                \     # exponential scaling factor in the Mahrt (1987) stability function
	)

# output: stability corrections
scalarRiBulkGround,               \  # bulk Richardson number for the ground surface (-)
scalarGroundStabilityCorrection,  \  # stability correction for the ground surface (-)
# output: scalar resistances
scalarFrictionVelocity,           \  # friction velocity (m s-1)
scalarGroundResistance,           \  # below canopy aerodynamic resistance (s m-1)
# output: derivatives in scalar resistances
dGroundResistance_dTGround,       \  # derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
							# Output aerodynamic resistances
							= aeroOut

# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# ***** TURBULENT HEAT FLUXES  **************************************************************************************************************************************
# *******************************************************************************************************************************************************************
# *******************************************************************************************************************************************************************

# check the need to compute numerical derivatives
if(ix_fDerivMeth == numerical)then
	nFlux=5  # compute the derivatives using one-sided finite differences
else:
	nFlux=1  # compute analytical derivatives

# either one or multiple flux calls, depending on if using analytical or numerical derivatives
for itry in np.arrange(nFlux,1,-1):  # (work backwards to ensure all computed fluxes come from the un-perturbed case)

	# -------------------------------------------------------------------------------------
	# state perturbations for numerical deriavtives with one-sided finite differences
	# note: no perturbations performed using analytical derivatives (nFlux=1)
	# -------------------------------------------------------------------------------------

	# identify the type of perturbation
	select case(itry)

	# un-perturbed case
	case(unperturbed)
	groundTemp        = groundTempTrial

	# perturb ground temperature
	case(perturbStateGround)
	groundTemp        = groundTempTrial + dx

	# perturb canopy air temperature
	groundTemp        = groundTempTrial

	# check for an unknown perturbation
	case default; err=10; message=trim(message)//"unknown perturbation"; return

	end select # (type of perturbation)

	# compute the saturation vapor pressure for ground temperature
	# NOTE: saturated vapor pressure derivatives don't seem that accurate....
	TG_celcius = groundTemp - Tfreeze
	call satVapPress(TG_celcius, scalarSatVP_GroundTemp, dSVPGround_dGroundTemp)

	# -------------------------------------------------------------------------------------
	# calculation block (unperturbed fluxes returned [computed last])
	# -------------------------------------------------------------------------------------

	# re-compute aerodynamic resistances for perturbed cases
	# NOTE: unperturbed fluxes computed earlier, and not over-written
	if itry not unperturbed:
		aeroOut = aeroResist(/
			# input: model control
			.false.,                                 / # intent(in): logical flag if would like to compute analytical derivaties
			ix_astability,                           / # intent(in): choice of stability function
			# input: above-canopy forcing data
			mHeight,                                 / # intent(in): measurement height (m)
			airtemp,                                 / # intent(in): air temperature at some height above the surface (K)
			windspd,                                 / # intent(in): wind speed at some height above the surface (m s-1)
			# input: temperature (canopy, ground, canopy air space)
			groundTemp,                              / # intent(in): ground temperature (K)
			# input: diagnostic variables
			scalarSnowDepth,                         / # intent(in): snow depth (m)
			# input: parameters
			z0Ground,                                / # intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
			critRichNumber,                          / # intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
			Louis79_bparam,                          / # intent(in): parameter in Louis (1979) stability function
			Mahrt87_eScale,                          / # intent(in): exponential scaling factor in the Mahrt (1987) stability function
		)
		
		# output: stability corrections
		notUsed_RiBulkGround,                    / # intent(out): bulk Richardson number for the ground surface (-)
		notUsed_scalarGroundStabilityCorrection, / # intent(out): stability correction for the ground surface (-)
		# output: scalar resistances
		notUsed_FrictionVelocity,                / # intent(out): friction velocity (m s-1)
		trialGroundResistance,                   / # intent(out): below canopy aerodynamic resistance (s m-1)
		# output: derivatives in scalar resistances
		notUsed_dGroundResistance_dTGround,      / # intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
									# Output aerodynamic resistance	
									= aeroOut

	# assign scalar resistances for un-perturbed cases
	else:
		trialGroundResistance = scalarGroundResistance

	# compute the relative humidity in the top soil layer and the resistance at the ground surface
	# NOTE: computations are based on start-of-step values, so only compute for the first flux call
	if firstFluxCall:
		# (soil water evaporation factor [0-1])
		soilEvapFactor = mLayerVolFracLiq(1)/(theta_sat - theta_res)
		# (resistance from the soil [s m-1])
		scalarSoilResistance = scalarGroundSnowFraction*1. + (1. - scalarGroundSnowFraction)*np.exp(8.2. - 4.22.*soilEvapFactor)  # Sellers (1992)
		# (relative humidity in the soil pores [0-1])
		if mLayerMatricHead(1) > -1*10**6: # avoid problems with numerical precision when soil is very dry
			soilRelHumidity_noSnow = exp( (mLayerMatricHead(1)*gravity) / (groundTemp*R_wv) )
		else
			soilRelHumidity_noSnow = 0.:
		scalarSoilRelHumidity  = scalarGroundSnowFraction*1. + (1. - scalarGroundSnowFraction)*soilRelHumidity_noSnow

	# compute turbulent heat fluxes
	turbOut turbFluxes(/
		# input: model control
		ix_fDerivMeth,                        / # intent(in): method used to calculate flux derivatives
		# input: above-canopy forcing data
		airtemp,                              / # intent(in): air temperature at some height above the surface (K)
		airpres,                              / # intent(in): air pressure of the air above the vegetation canopy (Pa)
		scalarVPair,                          / # intent(in): vapor pressure of the air above the vegetation canopy (Pa)
		# input: latent heat of sublimation/vaporization
		scalarLatHeatSubVapGround,            / # intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
		# input: canopy/ground temperature and saturated vapor pressure
		groundTemp,                           / # intent(in): ground temperature (K)
		scalarSatVP_GroundTemp,               / # intent(in): saturation vapor pressure at the temperature of the ground (Pa)
		dSVPGround_dGroundTemp,               / # intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
		# input: diagnostic variables
		scalarSoilRelHumidity,                / # intent(in): relative humidity in the soil pores [0-1]
		scalarSoilResistance,                 / # intent(in): resistance from the soil (s m-1)
		# input: derivatives in scalar resistances
		dGroundResistance_dTGround,           / # intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)

	# output: conductances (used to check derivative calculations)
	scalarGroundConductanceSH,            / # intent(out): ground conductance for sensible heat (m s-1)
	scalarGroundConductanceLH,            / # intent(out): ground conductance for latent heat -- includes soil resistance (m s-1)
	scalarEvapConductance,                / # intent(out): conductance for evaporation (m s-1)
	scalarTransConductance,               / # intent(out): conductance for transpiration (m s-1)
	scalarTotalConductanceSH,             / # intent(out): total conductance for sensible heat (m s-1)
	scalarTotalConductanceLH,             / # intent(out): total conductance for latent heat (m s-1)
	# output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
	scalarSenHeatGround,                  & # intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
	scalarLatHeatGround,                  & # intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
	# output: total heat fluxes to the atmosphere
	scalarSenHeatTotal,                   & # intent(out): total sensible heat flux to the atmosphere (W m-2)
	scalarLatHeatTotal,                   & # intent(out): total latent heat flux to the atmosphere (W m-2)
	# output: net fluxes
	turbFluxGround,                       & # intent(out): net turbulent heat fluxes at the ground surface (W m-2)
	# output: energy flux derivatives
	dTurbFluxGround_dTGround,             & # intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
	# output: liquid flux derivatives
	dLatHeatCanopyEvap_dTGround,          & # intent(out): derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)

	#notUsed_scalarCanopyStabilityCorrection  # stability correction for the canopy (-)
	#notUsed_scalarGroundStabilityCorrection  # stability correction for the ground surface (-)
	#notUsed_EddyDiffusCanopyTop              # eddy diffusivity for heat at the top of the canopy (m2 s-1)
	#notUsed_FrictionVelocity                 # friction velocity (m s-1)
	#notUsed_WindspdCanopyTop                 # windspeed at the top of the canopy (m s-1)
	#notUsed_WindspdCanopyBottom              # windspeed at the height of the bottom of the canopy (m s-1)
	#trialLeafResistance                      # mean leaf boundary layer resistance per unit leaf area (s m-1)
	#trialGroundResistance                    # below canopy aerodynamic resistance (s m-1)
	#trialCanopyResistance                    # above canopy aerodynamic resistance (s m-1)

	# save perturbed fluxes
	if ix_fDerivMeth == numerical:
		select case(itry) # (select type of perturbation)
		case(unperturbed)
		try0 = turbFluxGround
		exit
		
		case(perturbStateGround)
		try1 = turbFluxGround
		turbFluxGround_dStateGround = turbFluxGround          # total turbulent heat fluxes from the ground to the canopy air space (W m-2)
		latHeatCanEvap_dStateGround = scalarLatHeatCanopyEvap # perturbed value for the latent heat associated with canopy evaporation (W m-2)
		
		case default; err=10; message=trim(message)//"unknown perturbation"; return
		end select # (type of perturbation)

# test derivative
#if(ix_fDerivMeth == numerical)  print*, 'try0, try1 = ', try0, try1
#if(ix_fDerivMeth == numerical)  print*, 'derivative = ', (ix_fDerivMeth == numerical), (try1 - try0)/dx
#if(ix_fDerivMeth == analytical) print*, 'derivative = ', (ix_fDerivMeth == numerical), dTurbFluxGround_dTGround
#pause

# compute numerical derivatives
if(ix_fDerivMeth == numerical)then
	# derivatives w.r.t. ground temperature
	dTurbFluxGround_dTGround    = (turbFluxGround_dStateGround - turbFluxGround) / dx          # derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
	dLatHeatCanopyEvap_dTGround = (latHeatCanEvap_dStateGround - scalarLatHeatCanopyEvap) / dx # derivative in latent heat of canopy evaporation w.r.t. ground temperature (W m-2 K-1)
	dLatHeatCanopyEvap_dCanLiq = (latHeatCanEvap_dStateCanliq  - scalarLatHeatCanopyEvap) / dx # derivative in latent heat of canopy evaporation w.r.t. canopy liquid water content (J kg-1 s-1)

# test
#print*, (ix_fDerivMeth == numerical)
#print*, 'dTurbFluxCanair_dTCanair = ', dTurbFluxCanair_dTCanair
#print*, 'dTurbFluxCanair_dTCanopy = ', dTurbFluxCanair_dTCanopy
#print*, 'dTurbFluxCanair_dTGround = ', dTurbFluxCanair_dTGround
#print*, 'dTurbFluxCanopy_dTCanair = ', dTurbFluxCanopy_dTCanair
#print*, 'dTurbFluxCanopy_dTCanopy = ', dTurbFluxCanopy_dTCanopy
#print*, 'dTurbFluxCanopy_dTGround = ', dTurbFluxCanopy_dTGround
#print*, 'dTurbFluxGround_dTCanair = ', dTurbFluxGround_dTCanair
#print*, 'dTurbFluxGround_dTCanopy = ', dTurbFluxGround_dTCanopy
#print*, 'dTurbFluxGround_dTGround = ', dTurbFluxGround_dTGround
#print*, 'dLatHeatCanopyEvap_dCanLiq = ', dLatHeatCanopyEvap_dCanLiq
#print*, 'dLatHeatCanopyEvap_dTCanair = ', dLatHeatCanopyEvap_dTCanair
#print*, 'dLatHeatCanopyEvap_dTCanopy = ', dLatHeatCanopyEvap_dTCanopy
#print*, 'dLatHeatCanopyEvap_dTGround = ', dLatHeatCanopyEvap_dTGround
#print*, 'dTurbFluxCanair_dCanLiq = ', dTurbFluxCanair_dCanLiq
#print*, 'dTurbFluxCanopy_dCanLiq = ', dTurbFluxCanopy_dCanLiq
#print*, 'dTurbFluxGround_dCanLiq = ', dTurbFluxGround_dCanLiq
#pause

# (ground evaporation/sublimation)
if scalarLatHeatSubVapGround > LH_vap+verySmall: # sublimation
	# NOTE: this should only occur when we have formed snow layers, so check
	if nSnow == 0:
		error('only expect snow sublimation when we have formed some snow layers')
	scalarGroundEvaporation = 0.  # ground evaporation is zero once the snowpack has formed
	scalarSnowSublimation   = scalarLatHeatGround/LH_sub
	else
	# NOTE: this should only occur when we have no snow layers, so check
	if nSnow > 0 
		error('only expect ground evaporation when there are no snow layers')
	scalarGroundEvaporation = scalarLatHeatGround/LH_vap
	scalarSnowSublimation   = 0.  # no sublimation from snow if no snow layers have formed

# *******************************************************************************************************
# private subroutine aeroResist: compute aerodynamic resistances
# *******************************************************************************************************
def aeroResist(/
				# input: model control
				derivDesired,                  / # intent(in): flag to indicate if analytical derivatives are desired
				ixStability,                   / # intent(in): choice of stability function
				# input: above-canopy forcing data
				mHeight,                       / # intent(in): measurement height (m)
				airtemp,                       / # intent(in): air temperature at some height above the surface (K)
				windspd,                       / # intent(in): wind speed at some height above the surface (m s-1)
				# input: temperature (canopy, ground, canopy air space)
				groundTemp,                    / # intent(in): ground temperature (K)
				# input: diagnostic variables
				snowDepth,                     / # intent(in): snow depth (m)
				# input: parameters
				z0Ground,                      / # intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
				critRichNumber,                / # intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
				Louis79_bparam,                / # intent(in): parameter in Louis (1979) stability function
				Mahrt87_eScale,                / # intent(in): exponential scaling factor in the Mahrt (1987) stability function
				)
# -----------------------------------------------------------------------------------------------------------------------------------------
# compute aerodynamic resistances
# Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
#       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
#       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
# -----------------------------------------------------------------------------------------------------------------------------------------

	# local variables: general
	C_r = 0.3                     # roughness element drag coefficient (-) from Raupach (BLM, 1994)
	C_s = 0.00                # substrate surface drag coefficient (-) from Raupach (BLM, 1994)
	approxDragCoef_max = 0.   # maximum value of the approximate drag coefficient (-) from Raupach (BLM, 1994)

	# -----------------------------------------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------------------------------------
	# * compute resistance for the case without a canopy (bare ground, or canopy completely buried with snow)

	# check that measurement height above the ground surface is above the roughness length
	if mHeight < snowDepth+z0Ground:
		erro('measurement height < snow depth + roughness length')

	# compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)
	groundExNeut = (vkc**2.) / ( log((mHeight - snowDepth)/z0Ground)**2.) # turbulent transfer coefficient under conditions of neutral stability (-)
	groundResistanceNeutral = 1. / (groundExNeut*windspd)

	# define height above the snow surface
	heightAboveGround  = mHeight - snowDepth

	# check that measurement height above the ground surface is above the roughness length
	if heightAboveGround < z0Ground:
		print(\
				'z0Ground = %d'\n \
				'mHeight  = %d'\n \
				'snowDepth = %d'\n \
				'heightAboveGround = %d'\n \
				, (z0Ground,mHeight,snowDepth,heightAboveGround))
		error('height above ground < roughness length [likely due to snow accumulation]')

	# compute ground stability correction
	aStabilityOut = aStability(/
				   # input
				  derivDesired,                                     / # input: logical flag to compute analytical derivatives
				  ixStability,                                      / # input: choice of stability function
				  # input: forcing data, diagnostic and state variables
				  heightAboveGround,                                / # input: measurement height above the ground surface (m)
				  airtemp,                                          / # input: temperature above the ground surface (K)
				  groundTemp,                                       / # input: trial value of surface temperature -- "surface" is either canopy or ground (K)
				  windspd,                                          / # input: wind speed above the ground surface (m s-1)
				  # input: stability parameters
				  critRichNumber,                                   / # input: critical value for the bulk Richardson number where turbulence ceases (-)
				  Louis79_bparam,                                   / # input: parameter in Louis (1979) stability function
				  Mahrt87_eScale,                                   / # input: exponential scaling factor in the Mahrt (1987) stability function
				  )

	# output
	RiBulkGround,                                     / # output: bulk Richardson number (-)
	groundStabilityCorrection,                        / # output: stability correction for turbulent heat fluxes (-)
	dGroundStabilityCorrection_dRich,                 / # output: derivative in stability correction w.r.t. Richardson number for the ground surface (-)
	dGroundStabilityCorrection_dAirTemp,              / # output: (not used) derivative in stability correction w.r.t. air temperature (K-1)
	dGroundStabilityCorrection_dSfcTemp,              / # output: derivative in stability correction w.r.t. surface temperature (K-1)
										= aStabilityOut

	# compute the ground resistance (after stability corrections)
	groundResistance = groundResistanceNeutral/groundStabilityCorrection
	if groundResistance < 0.:
		error('ground resistance < 0 [no vegetation]'; return; endif

	# -----------------------------------------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------------------------------------
	# -----------------------------------------------------------------------------------------------------------------------------------------
	# ***** compute resistances for non-vegetated surfaces (e.g., snow)

	if derivDesired: # if analytical derivatives are desired

		# set canopy derivatives to zero (non-vegetated, remember)
		dCanopyResistance_dTCanopy = 0.
		dGroundResistance_dTCanopy = 0.

		# compute derivatives for ground resistance
		dGroundResistance_dTGround = -dGroundStabilityCorrection_dSfcTemp/(windspd*groundExNeut*groundStabilityCorrection**2.)

# * analytical derivatives not desired
# Set values to nan and return

# ********************************************************************************
# private subroutine turbFluxes: compute turbulent heat fluxes
# ********************************************************************************
def turbFluxes(&
				   # input: model control
				   computeVegFlux,                & # intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
				   ixDerivMethod,                 & # intent(in): choice of method used to compute derivative (analytical or numerical)
				   # input: above-canopy forcing data
				   airtemp,                       & # intent(in): air temperature at some height above the surface (K)
				   airpres,                       & # intent(in): air pressure of the air above the vegetation canopy (Pa)
				   VPair,                         & # intent(in): vapor pressure of the air above the vegetation canopy (Pa)
				   # input: latent heat of sublimation/vaporization
				   latHeatSubVapCanopy,           & # intent(in): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
				   latHeatSubVapGround,           & # intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
				   # input: canopy and ground temperature
				   canairTemp,                    & # intent(in): temperature of the canopy air space (K)
				   canopyTemp,                    & # intent(in): canopy temperature (K)
				   groundTemp,                    & # intent(in): ground temperature (K)
				   satVP_CanopyTemp,              & # intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
				   satVP_GroundTemp,              & # intent(in): saturation vapor pressure at the temperature of the ground (Pa)
				   dSVPCanopy_dCanopyTemp,        & # intent(in): derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
				   dSVPGround_dGroundTemp,        & # intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
				   # input: diagnostic variables
				   exposedVAI,                    & # intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
				   canopyWetFraction,             & # intent(in): fraction of canopy that is wet [0-1]
				   dCanopyWetFraction_dWat,       & # intent(in): derivative in the canopy wetted fraction w.r.t. total water content (kg-1 m-2)
				   dCanopyWetFraction_dT,         & # intent(in): derivative in wetted fraction w.r.t. canopy temperature (K-1)
				   canopySunlitLAI,               & # intent(in): sunlit leaf area (-)
				   canopyShadedLAI,               & # intent(in): shaded leaf area (-)
				   soilRelHumidity,               & # intent(in): relative humidity in the soil pores [0-1]
				   soilResistance,                & # intent(in): resistance from the soil (s m-1)
				   leafResistance,                & # intent(in): mean leaf boundary layer resistance per unit leaf area (s m-1)
				   groundResistance,              & # intent(in): below canopy aerodynamic resistance (s m-1)
				   canopyResistance,              & # intent(in): above canopy aerodynamic resistance (s m-1)
				   stomResistSunlit,              & # intent(in): stomatal resistance for sunlit leaves (s m-1)
				   stomResistShaded,              & # intent(in): stomatal resistance for shaded leaves (s m-1)
				   # input: derivatives in scalar resistances
				   dGroundResistance_dTGround,    & # intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
				   dGroundResistance_dTCanopy,    & # intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
				   dGroundResistance_dTCanair,    & # intent(in): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
				   dCanopyResistance_dTCanopy,    & # intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
				   dCanopyResistance_dTCanair,    & # intent(in): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
# -----------------------------------------------------------------------------------------------------------------------------------------
# local variables -- general
evapSmooth=1.             # smoothing parameter for latent heat (W m-2)
# local variables -- "constants"
# volHeatCapacityAir           # volumetric heat capacity of air (J m-3)
# latentHeatConstant           # latent heat constant (kg m-3 K-1)
# local variables -- derivatives for energy conductances
# dGroundCondSH_dGroundTemp    # derivative in ground conductance of sensible heat w.r.t. ground temperature
# local variables -- derivatives for mass conductances
# dGroundCondLH_dGroundTemp    # derivative in ground conductance w.r.t. ground temperature
# local variables -- sensible heat flux derivatives
# dSenHeatTotal_dTGround       # derivative in the total sensible heat flux w.r.t. ground temperature
# dSenHeatGround_dTGround      # derivative in the ground sensible heat flux w.r.t. ground temperature
# local variables -- latent heat flux derivatives
# dLatHeatGround_dTGround      # derivative in the ground latent heat flux w.r.t. ground temperature
# -----------------------------------------------------------------------------------------------------------------------------------------

# compute constants
volHeatCapacityAir = iden_air*cp_air           # volumetric heat capacity of air (J m-3)
latentHeatConstant = iden_air*w_ratio/airpres  # latent heat constant for (kg m-3 Pa-1)

# *****
# * compute conductances, and derivatives...
# ******************************************

# compute conductances for sensible heat (m s-1)
if computeVegFlux:
	leafConductance    = exposedVAI/leafResistance
	leafConductanceTr  = canopySunlitLAI/(leafResistance+stomResistSunlit) + canopyShadedLAI/(leafResistance+stomResistShaded)
	canopyConductance  = 1./canopyResistance
else:
	leafConductance    = 0.
	canopyConductance  = 0.
groundConductanceSH = 1./groundResistance

# compute total conductance for sensible heat
totalConductanceSH  = leafConductance + groundConductanceSH + canopyConductance

# compute conductances for latent heat (m s-1)
if computeVegFlux:
	evapConductance    = canopyWetFraction*leafConductance
	transConductance   = (1. - canopyWetFraction) * leafConductanceTr
else:
	evapConductance    = 0.
	transConductance   = 0.
groundConductanceLH = 1./(groundResistance + soilResistance)  # NOTE: soilResistance accounts for fractional snow, and =0 when snow cover is 100%
totalConductanceLH  = evapConductance + transConductance + groundConductanceLH + canopyConductance

# * compute derivatives
# NOTE: it may be more efficient to compute these derivatives when computing resistances
if ixDerivMethod == analytical:
	dGroundCondSH_dGroundTemp = -dGroundResistance_dTGround/groundResistance**2.         # derivative in ground conductance w.r.t. ground temperature
	dGroundCondLH_dGroundTemp = -dGroundResistance_dTGround/(groundResistance+soilResistance)**2. # derivative in ground conductance w.r.t. ground temperature

# *****
# * compute sensible and latent heat fluxes, and derivatives...
# *************************************************************

# compute sensible and latent heat fluxes from the ground to the canopy air space (W m-2)
senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - airtemp)                                                 # (positive downwards)
latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VPair)         # (positive downwards)
senHeatTotal       = senHeatGround

# compute latent heat flux from the canopy air space to the atmosphere
latHeatTotal = latHeatGround

# * compute derivatives
if ixDerivMethod == analytical:
# compute derivatives for the ground fluxes w.r.t. ground temperature
	dSenHeatGround_dTGround = (-volHeatCapacityAir*dGroundCondSH_dGroundTemp)*(groundTemp - airtemp) + &                                               # d(ground sensible heat flux)/d(ground temp)
							 (-volHeatCapacityAir*groundConductanceSH)
	dLatHeatGround_dTGround = (-latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp)*(satVP_GroundTemp*soilRelHumidity - VPair) + &       # d(ground latent heat flux)/d(ground temp)
							 (-latHeatSubVapGround*latentHeatConstant*groundConductanceLH)*dSVPGround_dGroundTemp*soilRelHumidity


# *****
# * compute net turbulent fluxes, and derivatives...
# **************************************************

# compute net fluxes
turbFluxCanair = senHeatTotal - senHeatCanopy - senHeatGround            # net turbulent flux at the canopy air space (W m-2)
turbFluxCanopy = senHeatCanopy + latHeatCanopyEvap + latHeatCanopyTrans  # net turbulent flux at the canopy (W m-2)
turbFluxGround = senHeatGround + latHeatGround                           # net turbulent flux at the ground surface (W m-2)

# * compute derivatives
if ixDerivMethod == analytical:
	# (energy derivatives)
	dTurbFluxGround_dTGround = dSenHeatGround_dTGround + dLatHeatGround_dTGround                                     # derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
else # (just make sure we return something)
	# (energy derivatives)
	dTurbFluxGround_dTGround = 0.


# *******************************************************************************************************
# private subroutine aStability: compute stability corrections for turbulent heat fluxes (-)
# *******************************************************************************************************
def aStability(/
				   # input: control
				   computeDerivative,              / # input: logical flag to compute analytical derivatives
				   ixStability,                    / # input: choice of stability function
				   # input: forcing data, diagnostic and state variables
				   mHeight,                        / # input: measurement height (m)
				   airTemp,                        / # input: air temperature (K)
				   sfcTemp,                        / # input: surface temperature (K)
				   windspd,                        / # input: wind speed (m s-1)
				   # input: stability parameters
				   critRichNumber,                 / # input: critical value for the bulk Richardson number where turbulence ceases (-)
				   Louis79_bparam,                 / # input: parameter in Louis (1979) stability function
				   Mahrt87_eScale,                 / # input: exponential scaling factor in the Mahrt (1987) stability function
				   # output
				   RiBulk,                         / # output: bulk Richardson number (-)
				   stabilityCorrection,            / # output: stability correction for turbulent heat fluxes (-)
				   dStabilityCorrection_dRich,     / # output: derivative in stability correction w.r.t. Richardson number (-)
				   dStabilityCorrection_dAirTemp,  / # output: derivative in stability correction w.r.t. temperature (K-1)
				   dStabilityCorrection_dSfcTemp,  / # output: derivative in stability correction w.r.t. temperature (K-1)
				   )
# compute the bulk Richardson number (-)
bulkRichardsonOut = bulkRichardson(/
								# input
								airTemp,                        / # input: air temperature (K)
								sfcTemp,                        / # input: surface temperature (K)
								windspd,                        / # input: wind speed (m s-1)
								mHeight,                        / # input: measurement height (m)
								computeDerivative,              / # input: flag to compute the derivative
								)

# output
RiBulk,                         / # output: bulk Richardson number (-)
dRiBulk_dAirTemp,               / # output: derivative in the bulk Richardson number w.r.t. air temperature (K-1)
dRiBulk_dSfcTemp,               / # output: derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
							= bulkRichardonOut
							
# set derivative to one if not computing it
if not computeDerivative:
	dStabilityCorrection_dRich    = 1.
	dStabilityCorrection_dAirTemp = 1.
	dStabilityCorrection_dSfcTemp = 1.

# ***** process unstable cases
if RiBulk<0.:
	# compute surface-atmosphere exchange coefficient (-)
	stabilityCorrection = (1. - 16.*RiBulk)**0..
# compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
if computeDerivative:
	dStabilityCorrection_dRich    = (-16.) * 0..*(1. - 16.*RiBulk)**(-0..)
	dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
	dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich

# ***** process stable cases
select case(ixStability)

#### ("standard" stability correction, a la Anderson 1976)
case(standard)
# compute surface-atmosphere exchange coefficient (-)
if RiBulk <  critRichNumber:
	stabilityCorrection = (1. - 5.*RiBulk)**2.
if RiBulk >= critRichNumber:
	####### EPSILON IS NOT A DEFINED FUNCTION
	stabilityCorrection = epsilon(stabilityCorrection)
# compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
if computeDerivative:
	if RiBulk <  critRichNumber:
		dStabilityCorrection_dRich = (-5.) * 2.*(1. - 5.*RiBulk)
	if RiBulk >= critRichNumber:
		dStabilityCorrection_dRich = 0.

# (Louis 1979)
case(louisInversePower)
# scale the "b" parameter for stable conditions
bprime = Louis79_bparam/2.
# compute surface-atmosphere exchange coefficient (-)
stabilityCorrection = 1. / ( (1. + bprime*RiBulk)**2. )
if stabilityCorrection < epsilon(stabilityCorrection):
	stabilityCorrection = epsilon(stabilityCorrection)
# compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
if computeDerivative:
	dStabilityCorrection_dRich = bprime * (-2.)*(1. + bprime*RiBulk)**(-3.)

# (Mahrt 1987)
case(mahrtExponential)
# compute surface-atmosphere exchange coefficient (-)
stabilityCorrection = exp(-Mahrt87_eScale * RiBulk)
if stabilityCorrection < epsilon(stabilityCorrection):
	stabilityCorrection = epsilon(stabilityCorrection)
# compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
if computeDerivative:
	dStabilityCorrection_dRich = (-Mahrt87_eScale) * exp(-Mahrt87_eScale * RiBulk)

# (return error if the stability correction method is not found)
case default
err=10; message=trim(message)//"optionNotFound[stability correction]"; return

# endselect

# get the stability correction with respect to air temperature and surface temperature
# NOTE: air temperature is used for canopy air temperature, which is a model state variable
if computeDerivative:
	dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
	dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich

# *******************************************************************************************************
# private subroutine bulkRichardson: compute bulk Richardson number
# *******************************************************************************************************
def bulkRichardson(\
				   # input
				   airTemp,                    \ # input: air temperature (K)
				   sfcTemp,                    \ # input: surface temperature (K)
				   windspd,                    \ # input: wind speed (m s-1)
				   mHeight,                    \ # input: measurement height (m)
				   computeDerivative,          \ # input: flag to compute the derivative
					)
# compute local variables
T_grad = airtemp - sfcTemp
T_mean = 0..*(airtemp + sfcTemp)
RiMult = (gravity*mHeight)/(windspd*windspd)
# compute the Richardson number
RiBulk = (T_grad/T_mean) * RiMult
# compute the derivative in the Richardson number
if computeDerivative:
	dRiBulk_dAirTemp =  RiMult/T_mean - RiMult*T_grad/(0..*((airtemp + sfcTemp)**2.))
	dRiBulk_dSfcTemp = -RiMult/T_mean - RiMult*T_grad/(0..*((airtemp + sfcTemp)**2.))
else:
	dRiBulk_dAirTemp = 1.
	dRiBulk_dSfcTemp = 1.


