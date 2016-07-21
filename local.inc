' This file is included from TOPAS.INC
' Place local standard macros here

macro Out_cif_powder_data(file)
{
   xdd_out file load out_record out_fmt out_eqn
   {
		 " %11.5f  " = X;
       " %11.5f  " = Yobs;
       " %11.5f  " = Ycalc;
       " %11.5f\n" = 1/Yobs;
   }
}

'Description: ride H atoms in organic structure
'Contributed by: Alan Coelho
macro H_ride(h_site, c_site, xv, yv, zv, beqcv)
{
	site h_site
		x = Constant((xv)-Get(Find_Child(Get(site_recs), c_site), x)) + Get(Find_Child(Get(site_recs), c_site), x); : xv
		y = Constant((yv)-Get(Find_Child(Get(site_recs), c_site), y)) + Get(Find_Child(Get(site_recs), c_site), y); : yv
		z = Constant((zv)-Get(Find_Child(Get(site_recs), c_site), z)) + Get(Find_Child(Get(site_recs), c_site), z); : zv
		occ H 1
		beq beqcv
}



macro Distance_Restrain_Morse(sites, t, t_calc, a, wscale)
{
	Distance_Restrain_Morse(sites,, t, t_calc, a, wscale)
}

macro Distance_Restrain_Morse(sites, c, t, t_calc, a, wscale)
{
	#m_ifarg c ""
      #m_unique_not_refine c
   #m_endif
	sites_distance c 
		load site_to_restrain { sites }
	prm = c; : t_calc
	'prm = wscale*(1-Exp(-a(c-t)^2)); : pen
    penalty = (wscale) (1-Exp(-a(c-t)))^2; 
}

macro Append_2Th_Ip_file(file, wI)
{
   phase_out file append load out_record out_fmt out_eqn
   {
      " %11.5f"   = 2 Rad Th;
      " %11.5f\n" = wI;
   }
}

macro Append_d_Ip_file(file, wI)
{
   phase_out file append load out_record out_fmt out_eqn
   {
      " %11.5f"   = Lam / (2 Sin(Th));
      " %11.5f\n" = wI;
   }
}


macro Append_2Th_Ip_file(file) {   Append_2Th_Ip_file(file, I_no_scale_pks) }
macro Append_d_Ip_file(file) { Append_d_Ip_file(file, I_no_scale_pks) }
macro Append_2Th_IScaled_file(file) {   Append_2Th_Ip_file(file, I_after_scale_pks) }
macro Append_d_IScaled_file(file) { Append_d_Ip_file(file, I_after_scale_pks) }


		
		macro Out_Graphs(file)
{
   xdd_out file load out_record out_fmt out_eqn
   {
	" %11.5f  " = X;       
	" %11.5f  " = Yobs;
       " %11.5f  " = Ycalc;
       " %11.5f\n" = Yobs-Ycalc;
	     }
}

macro Distance_Restrain_Breakable(sites, t, t_calc, a, wscale)
{
	Distance_Restrain_Breakable(sites,, t, t_calc, a, wscale)
}



macro Distance_Restrain_Breakable(sites, c, t, t_calc, a, wscale)
{
	#m_ifarg c ""
      #m_unique_not_refine c
   #m_endif
	sites_distance c 
		load site_to_restrain { sites }
	prm = c; : t_calc
	'prm = wscale*(1-Exp(-a(c-t)^2)); : pen
    penalty = (wscale) (1-Exp(-a(Abs(c-t))))^2; 
}




macro Stephens_tetragonal(s400, s004, s220, s202, eta)
{
prm mhkl = Abs(s400 (H^4+ K^4)+ s004 L^4+ s220 H^2 K^2 + s202 (H^2 L^2 + K^2 L^2) );
 
prm pp = D_spacing^2 * Sqrt(Max(mhkl,0)) / 1000;
gauss_fwhm = 1.8/3.1415927 pp (1-eta) Tan(Th) + 0.0001;
lor_fwhm = 1.8/3.1415927 pp eta Tan(Th) + 0.0001;
}


macro Stephens_monoclinic(s400, s040, s004, s220, s202, s022, s301, s121, s103, eta)
{
prm mhkl = H^4 s400 + K^4 s040 + L^4 s004 +
H^2 K^2 s220 + H^2 L^2 s202 + K^2 L^2 s022 +
H K^2 L s121 +
H L^3 s103 + H^3 L s301;
 
prm pp = D_spacing^2 * Sqrt(Max(mhkl,0)) / 1000;
 
gauss_fwhm = 1.8/3.1415927 pp (1-eta) Tan(Th) + 0.0001;
lor_fwhm = 1.8/3.1415927 pp eta Tan(Th) + 0.0001;
}

macro Stephens_hexagonal(s400, s202, s004, eta)
{
prm mhkl = H^4 s400 + K^4 s400 + L^4 s004 +
H^2 K^2 3 s400 + H^2 L^2 s202 + K^2 L^2 s202 +
H K L^2 s202 +
H^3 K 2 s400 + H K^3 2 s400;
 
prm pp = D_spacing^2 * Sqrt(Max(mhkl,0)) / 1000;
 
gauss_fwhm = 1.8/3.1415927 pp (1-eta) Tan(Th) + 0.0001;
lor_fwhm = 1.8/3.1415927 pp eta Tan(Th) + 0.0001;
}



macro Stephens_orthorhombic(s400, s040, s004, s220, s202, s022, eta)
{
prm mhkl = H^4 s400 + K^4 s040 + L^4 s004 +
H^2 K^2 s220 + H^2 L^2 s202 + K^2 L^2 s022;
 
prm pp = D_spacing^2 * Sqrt(Max(mhkl,0)) / 1000;
 
gauss_fwhm = 1.8/3.1415927 pp (1-eta) Tan(Th) + 0.0001;
lor_fwhm = 1.8/3.1415927 pp eta Tan(Th) + 0.0001;
}


macro Foil_Thickness( v ) {Foil_Thickness(,v)}
macro Foil_Thickness(c, v)
{
   #m_argu c
   If_Prm_Eqn_Rpt(c, v, min 0 max 2)
   gauss_fwhm = 180/Pi CeV(c, v)/Rs Tan(Th);
}

macro Foil_ThicknessI( v ) {Foil_Thickness(,v)}
macro Foil_ThicknessI(c, v)
{
   #m_argu c
   If_Prm_Eqn_Rpt(c, v, min 0 max 2)
   gauss_fwhm = 180/Pi CeV(c, v)/Rs Tan(Th);
   scale_pks = 1/Cos(Th);
}


macro Angle_Restrain_Breakable(sites, t, t_calc, a, wscale)
{
	Angle_Restrain_Breakable(sites,, t, t_calc, a, wscale)
}



macro Angle_Restrain_Breakable(sites, c, t, t_calc, a, wscale)
{
	#m_ifarg c ""
      #m_unique_not_refine c
   #m_endif
	sites_angle c 
		load site_to_restrain { sites }
	prm = c; : t_calc
	'prm = wscale*(1-Exp(-a(c-t)^2)); : pen
    penalty = (wscale) (1-Exp(-a(Abs(c-t))))^2; 
}



macro Foil_ThicknessAbsI( thick, abs ) {Foil_ThicknessAbsI(,thick,,abs)}
macro Foil_ThicknessAbsI(ct, thick, ca, abs)
{
   #m_argu ct
   If_Prm_Eqn_Rpt(ct, thick, min 0 max 2)
   #m_argu ca
   If_Prm_Eqn_Rpt(ca, abs, min 1 max 1000)   
   
   hat = 180/Pi CeV(ct, thick)/Rs Tan(Th);
   ' scale_pks = 1/Cos(Th);
   scale_pks = 1/Cos(Th) Exp(-0.1 CeV(ca, abs) CeV(ct, thick)/Cos(Th))/Exp(-0.1 CeV(ca, abs) CeV(ct, thick));
}

macro Foil_ThicknessAbsIgauss( thick, abs ) {Foil_ThicknessAbsI(,thick,,abs)}
macro Foil_ThicknessAbsIgauss(ct, thick, ca, abs)
{
   #m_argu ct
   If_Prm_Eqn_Rpt(ct, thick, min 0 max 2)
   #m_argu ca
   If_Prm_Eqn_Rpt(ca, abs, min 1 max 1000)   
   
   gauss_fwhm = 180/Pi CeV(ct, thick)/Rs Tan(Th);
   ' scale_pks = 1/Cos(Th);
   scale_pks = 1/Cos(Th) Exp(-0.1 CeV(ca, abs) CeV(ct, thick)/Cos(Th))/Exp(-0.1 CeV(ca, abs) CeV(ct, thick));
}


macro Out_CIF_IUCR(namme) {
	macro Outname() { namme }
	Out_CIF_IUCR(Outname().cif, namme)
}

macro Out_CIF_IUCR(file, myname)
{
   out file
	  Out_String("\ndata_")
	  Out_String(myname)
      Out_String("\n_chemical_name_mineral ")
	  Out_String(myname)
	  Out_String("\n")
	  Out_String("#Copy this section from PLATON")
	  Out_String("
#===============================================================================

# 6. CRYSTAL DATA (replace F000 and xyz with PLATON data, fill exp. details)
")
      Out(Get(a), "\n_cell_length_a  %V")
      Out(Get(b), "\n_cell_length_b  %V")
      Out(Get(c), "\n_cell_length_c  %V")
      Out(Get(al), "\n_cell_angle_alpha %V")
      Out(Get(be), "\n_cell_angle_beta  %V")
      Out(Get(ga), "\n_cell_angle_gamma %V")
      Out(Get(cell_volume), "\n_cell_volume %V")
	  Out(Lam, "\n_diffrn_radiation_wavelength %V\n")
      Out(Get(sp_grp_char), "\n#Copy full symmetry info from PLATON\n_symmetry_space_group_name_H-M %s")
load out_record out_fmt out_eqn {
"\n_exptl_crystal_density_diffrn   %.3f" =Phase_Density_Eqn_g_on_cm3;
"\n_exptl_absorpt_coefficient_mu %.3f" =Phase_LAC_Eqn_1_on_cm/10;
}
	Out_String("\n_exptl_crystal_F_000     ? #PLATON")
      Out_String("\n\n       #REPLACE with PLATON (i.e. canonical) output to make tables correct\nloop_\n_symmetry_equiv_pos_as_xyz")
      Out(Get(sp_xyzs_txt),  "%s")
	  
	  Out_String("

 _cell_measurement_temperature 298

_pd_spec_mounting                          
; 
deposited on sticky tape 
;
_pd_spec_mount_mode         transmission    # options are reflection
                                          # or transmission

_pd_spec_shape              flat_sheet    # options are cylinder
                                          # flat_sheet OR irregular

_pd_char_particle_morphology      round
_pd_char_colour                   white   # use ICDD colour descriptions




_pd_prep_cool_rate                ?        
_pd_prep_pressure                 ?       
_pd_prep_temperature              ?
 
 
_exptl_absorpt_correction_type               analytical

_exptl_absorpt_process_details                
;
International Tables Vol C Tables 2.3.1.20
 A(Th) = (t/cos(Th)) exp(-s/cos(Th))
refined coefficients are
 s =  ??
 t =  ??
;
")	  

Out_String("#===============================================================================

# 7. EXPERIMENTAL DATA (leave as is)


_exptl_special_details                   ?

_diffrn_radiation_probe              x-ray
_diffrn_radiation_type               CuK\a~1~
_diffrn_radiation_source             ?
_diffrn_radiation_monochromator      
;
Ge 111
;
_diffrn_ambient_temperature     298                    
_diffrn_source              
;
sealed X-ray tube
;        
_diffrn_source_target            Cu                                     
_diffrn_measurement_device_type  
;
Bruker D8 Advance Vario
;       
_diffrn_detector               
;
1D linear
;
_diffrn_detector_type          
;
Bruker LynxEye compound silicon strip 1D linear detector 
;
_pd_meas_scan_method              cont        

_computing_data_collection           
;
Bruker XRD Commander 2.6 (Bruker, 2007)
;
_computing_cell_refinement           
;
TOPAS 4.2 (Coelho, 2003; Bruker, 2009)
;
_computing_data_reduction            
;
TOPAS 4.2 (Coelho, 2003; Bruker, 2009)
;
_computing_structure_solution        
;
FOX 1.9.7.0 (Favre-Nicolin, Cerny, 2002)
;
_computing_structure_refinement      
;
TOPAS 4.2 (Coelho, 2003; Bruker, 2009)
;
_computing_molecular_graphics        
;
OLEX2 1.2 (Dolomanov et al., 2009, 2)
;
_computing_publication_material      
;
PLATON (Spek, 2003)
;
")
Out_String("\n#Beware: these are values as defined by start_X and finish_X.\n# Real measurerement range may differ")
load out_record out_fmt out_eqn {
"\n_pd_meas_number_of_points %.0f" =(X2-X1)/Peak_Calculation_Step+1;
"\n_pd_meas_2theta_range_min  %.8f" =X1;
"\n_pd_meas_2theta_range_max %.8f" =X2;
"\n_pd_meas_2theta_range_inc %.8f" = Peak_Calculation_Step;
}
	  
	  Out_String("
#===============================================================================

# 8. REFINEMENT DATA (fill in the blanks, here PLATON will not help)
")
	  Out_String("\n_pd_proc_ls_special_details              ?")
		Out_String("\n_pd_proc_ls_profile_function      pseudo-Voigt")
		Out_String("\n#EDIT!\n_pd_proc_ls_background_function\n;\n??th order Chebyshev polynomial \n;")
Out_String("
_pd_proc_ls_pref_orient_corr     
;
EDIT!
?? order hkl dependent symmetrized spherical harmonics series (Jarvinen 1993).

The coefficients used were:

y00    1.000
y20 ...

;
")
	load out_record out_fmt out_eqn {
		"\n_pd_proc_ls_prof_R_factor         %.5f"  =Get(r_p)/100; 
		"\n_pd_proc_ls_prof_wR_factor   %.5f" =Get(r_wp)/100;
		"\n_pd_proc_ls_prof_wR_expected   %.5f" =Get(r_exp)/100;
		"\n_refine_ls_goodness_of_fit_all  %.3f" =Get(gof);
	}
	Out_String("
_refine_ls_R_I_factor ? #EDIT! This is r_bragg/100
")

	
	Out_String("
_refine_ls_structure_factor_coef     Inet
_refine_ls_matrix_type               fullcycle
_refine_ls_weighting_scheme          calc
_refine_ls_weighting_details         1/Y~i~
_refine_ls_hydrogen_treatment        constr
_refine_ls_extinction_method         none        
_refine_ls_extinction_coef           ?  
")
	Out_String("
_refine_ls_number_parameters         ? #EDIT!
_refine_ls_number_restraints         ? #EDIT!
_refine_ls_number_constraints         ? #EDIT!
")




' # The following four items apply to angular dispersive measurements.
' # 2theta minimum, maximum and increment (in degrees) are for the 
' # intensities used in the refinement. 


Out(X1, "\n_pd_proc_2theta_range_min %.6f")  
Out(X2, "\n_pd_proc_2theta_range_max %.6f")
Out(Peak_Calculation_Step, "\n_pd_proc_2theta_range_inc %.6f")


 

	  Out(Lam, "\n_pd_proc_wavelength %V\n")
		Out_String("\n_pd_block_diffractogram_id     ? #EDIT!")



		Out_String("\n\n#Do not forget to run fix_atom_prec.pl <outfile> <ciffile> \n# and replace everything below with the rest of PLATON file\nloop_")
		Out_String("\n_atom_site_label")
         Out_String("\n_atom_site_type_symbol")
		 Out_String("\n_atom_site_adp_type")
         Out_String("\n_atom_site_site_symmetry_multiplicity")
         Out_String("\n_atom_site_fract_x")
         Out_String("\n_atom_site_fract_y")
         Out_String("\n_atom_site_fract_z")
         Out_String("\n_atom_site_occupancy")
         Out_String("\n_atom_site_B_iso_or_equiv")
		 ' Out_String("\n_atom_site_U_iso_or_equiv")
         atom_out file append
            load out_record out_fmt out_eqn
            {
               "\n%-5s" = Get_From_String(Get(current_atom), site);
				" %-3s" = Get_From_String(Get(current_atom), atom);
			}	
			Out_String(" Biso ")
			load out_record out_fmt out_eqn
			{
               " %3.0f" = Get_From_String(Get(current_atom), num_posns);
				" %V" = Get_From_String(Get(current_atom), x);
               " %V" = Get_From_String(Get(current_atom), y);
               " %V" = Get_From_String(Get(current_atom), z);
               " %V" = Get_From_String(Get(current_atom), occ);
               " %V" = Get_From_String(Get(current_atom), beq);
            }
			' out_record
				' out_eqn = ((Get_From_String(Get(current_atom), beq))/(8*Pi^2));
				' out_fmt " %2.5f"
			' out_record
				' out_eqn = (Get_Prm_Error((Get_From_String(Get(current_atom), beq)))/(8*Pi^2));
				' out_fmt " %2.5f"
			
}

macro Toraya_Weighting(toraya_e) {
	weighting = 1 / (SigmaYobs^2*Yobs^(toraya_e-1));
}

macro Out_Graphs_Weight(file)
{
   xdd_out file load out_record out_fmt out_eqn
   {
	" %11.9f  " = X;       
	" %11.9f  " = Yobs;
       " %11.9f  " = Ycalc;
       " %11.9f  " = Yobs-Ycalc;
	   " %11.9f " = Get(weighting);
	   " %11.9f\n" = SigmaYobs;
	     }
}