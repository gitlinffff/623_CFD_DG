#!/bin/bash

# --- Configuration ---
EXE_PATH="/home/linfel/umich_course/623_DG/dg_crv/build_unstd/main"
MESH_DIR="/home/linfel/umich_course/623_DG/dg_crv/mesh"
BASE_DIR="/home/linfel/umich_course/623_DG/dg_crv_exp"
MESHES=("mesh2")
ORDERS=(0 1 2)

# --- CFL Logic ---
# Edit these values to set the CFL for each polynomial order
get_cfl() {
    local p_order=$1
    case $p_order in
        0) echo "0.6"   ;; # P0 is very stable
        1) echo "0.5"   ;; # P1
        2) echo "0.3"   ;; # P2 (Curved)
        3) echo "0.05"  ;; # P3 (Curved/High Order)
        *) echo "0.05"  ;; # Default
    esac
}

# --- Main Loop ---
for P in "${ORDERS[@]}"; do
    for M in "${MESHES[@]}"; do
        
        # Determine specific parameters for this case
        CURRENT_CFL=$(get_cfl $P)
        CASE_NAME="${M}_p${P}"
        CASE_DIR="${BASE_DIR}/${CASE_NAME}"
        
        echo "Setting up ${CASE_NAME} with CFL=${CURRENT_CFL}..."
        
        # 1. Create directory
        mkdir -p "${CASE_DIR}"
        
        # 2. Generate input.dat
        sed -e "s/{{MESH}}/${M}/g" \
            -e "s/{{ORDER}}/${P}/g" \
            -e "s/{{CFL}}/${CURRENT_CFL}/g" \
            "${BASE_DIR}/input.template" > "${CASE_DIR}/input.dat"
            
        # 3. Generate submit.sub
        sed -e "s/{{MESH}}/${M}/g" \
            -e "s/{{ORDER}}/${P}/g" \
            "${BASE_DIR}/submit.template" > "${CASE_DIR}/submit.sub"
            
        # 4. Create soft link
        ln -sf "${EXE_PATH}" "${CASE_DIR}/main"
				ln -sf "${MESH_DIR}" "${CASE_DIR}"
        
        # 5. Submit (Optional: uncomment line below to auto-submit)
        cd "${CASE_DIR}" && sbatch submit.sub && cd "${BASE_DIR}"
        
    done
done

echo "Done. Folders created in ${BASE_DIR}"
