#!/usr/bin/env bash

#chgrp adm -R .
#chmod g+w+x -R .
# Generate and build reactions

cd tcg_slb_database
./scripts/generate_reactions_eclogite -v slb21
./scripts/build_reactions database/reactions/eclogitization_2024_slb21_rx.rxml