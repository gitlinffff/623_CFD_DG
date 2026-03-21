#!/bin/bash

BASE_DIR=$(pwd)
MESHES=("mesh0" "mesh1" "mesh2")
ORDERS=(0 1)

for M in "${MESHES[@]}"; do
    for P in "${ORDERS[@]}"; do
        
        # Path to the results directory
        DIR_PATH="${M}_p${P}/results/${M}_curved_p${P}"
        
        if [ -d "$DIR_PATH/unsteady" ]; then
            SHORT_M=$(echo $M | sed 's/mesh/m/')
            FOLDER_NAME="${SHORT_M}p${P}_crv_unstd"
            ZIP_NAME="${FOLDER_NAME}.zip"
            
            echo "Zipping $FOLDER_NAME..."

            # 1. Enter the results directory
            cd "$DIR_PATH" || exit

            # 2. Rename 'unsteady' to the zip name temporarily
            mv unsteady "$FOLDER_NAME"

            # 3. Zip the renamed folder (creates the top-level folder inside zip)
            zip -rq "$BASE_DIR/$ZIP_NAME" "$FOLDER_NAME"

            # 4. Rename it back to 'unsteady'
            mv "$FOLDER_NAME" unsteady

            # 5. Go back to the root directory for the next loop
            cd "$BASE_DIR" || exit
            
        else
            echo "Skip: $M p$P results not found."
        fi
    done
done

echo "Finished! All zips are in $(pwd)"
