
if(PRESSIO_ENABLE_TPL_EIGEN)
########################
# hessian-gradient API
########################
add_serial_exe_and_test(gn_neq_frozen_jac_with_hessian_grad_api_eigen ${ROOTNAME} frozen_jac_with_hessian_grad_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_neq_frozen_jac_with_hessian_grad_api_eigen PRIVATE USE_GN_NEQ)

add_serial_exe_and_test(lm_neq_frozen_jac_with_hessian_grad_api_eigen ${ROOTNAME} frozen_jac_with_hessian_grad_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_lm_neq_frozen_jac_with_hessian_grad_api_eigen PRIVATE USE_LM_NEQ)

#############################
# fused hessian-gradient API
#############################
add_serial_exe_and_test(gn_neq_frozen_jac_with_fused_hessian_grad_api_eigen ${ROOTNAME} frozen_jac_with_fused_hessian_grad_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_neq_frozen_jac_with_fused_hessian_grad_api_eigen PRIVATE USE_GN_NEQ)

add_serial_exe_and_test(lm_neq_frozen_jac_with_fused_hessian_grad_api_eigen ${ROOTNAME} frozen_jac_with_fused_hessian_grad_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_lm_neq_frozen_jac_with_fused_hessian_grad_api_eigen PRIVATE USE_LM_NEQ)


########################
# residual-jacobian API
########################
add_serial_exe_and_test(gn_neq_frozen_jac_with_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_neq_frozen_jac_with_res_jac_api_eigen PRIVATE USE_GN_NEQ)

add_serial_exe_and_test(gn_qr_frozen_jac_with_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_qr_frozen_jac_with_res_jac_api_eigen PRIVATE USE_GN_QR)

add_serial_exe_and_test(lm_neq_frozen_jac_with_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_lm_neq_frozen_jac_with_res_jac_api_eigen PRIVATE USE_LM_NEQ)

##############################
# fused residual-jacobian API
##############################
add_serial_exe_and_test(gn_neq_frozen_jac_with_fused_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_fused_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_neq_frozen_jac_with_fused_res_jac_api_eigen PRIVATE USE_GN_NEQ)

add_serial_exe_and_test(gn_qr_frozen_jac_with_fused_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_fused_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_gn_qr_frozen_jac_with_fused_res_jac_api_eigen PRIVATE USE_GN_QR)

add_serial_exe_and_test(lm_neq_frozen_jac_with_fused_res_jac_api_eigen ${ROOTNAME} frozen_jac_with_fused_res_jac_api.cc "PASSED")
target_compile_definitions(${ROOTNAME}_lm_neq_frozen_jac_with_fused_res_jac_api_eigen PRIVATE USE_LM_NEQ)

endif()
