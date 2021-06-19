# Overview

This repository contains scripts to reproduce figures found in the [FAUST Manuscript](https://www.biorxiv.org/content/10.1101/702118v2).

To reproduce a given figure, required datasets must be placed in the cooresponding directories under `/dataSets/publication_gating_sets`.

Then, corresponding scripts under `/faustRuns` must be run.

Finally, scripts under `createFigures` must be run to re-create the figures.

Required R libraries are implied by the `library` statements at a head of a given script.

# Structure of `/dataSets/publication_gating_sets`.

This directory structure expected by the scripts is provided in the following diagram.

```
dataSets/publication_gating_sets/
|-- CITN-07
|   |-- citn07_longitudinal_gs
|-- CITN-09
|   |-- cryopreserved_pbmcs
|   |   |-- citn09_cryo_pbmcs_gs
|   |-- myeloid_panel
|   |   |-- citn09_myeloid_gs
|   `-- whole_blood_samples
|       |-- citn09_fresh_whole_blood_gs
|-- melanoma_cytof
|   |-- FlowRepository_FR-FCM-ZYKP_files
|   |   `-- attachments
|   `-- melanoma_gs
`-- melanoma_facs
    |-- melanoma_facs_gs
```

The named gating set data structures must be downloaded and placed into the location implied by the above diagram.

If a named gating set is not provided for download, it must be downloaded from [FlowRepository](https://flowrepository.org/id/FR-FCM-ZYKP).



