# GENIE Tunes and Model Configurations

## Overview

GENIE tunes define a consistent set of physics models and parameters used in neutrino interaction simulations. A tune is built on top of a **Comprehensive Model Configuration (CMC)**, which specifies the collection of physics models and their configurations. A tune additionally defines the parameters, fitting datasets, dataset weights, and prior constraints used to tune the configuration.

## Building a GENIE Tune

A GENIE tune consists of three main components:

1. **Model configuration (CMC):** Defines the physics models used in the simulation and their associated configurations (ModelConfiguration.xml).
2. **Parameter set:** Defines the parameters of interest and nuisance parameters to be tuned (CommonPram.xml).
Thie subfolders can contain modifications to those parameters.




## Model Consistency and Parameter Naming

**All models within a tune must use consistent configurations and parameter definitions.**

GENIE parameters are identified by their names. Consequently, parameter names must be unique and used consistently across all models. If different models use the same parameter name for different physical quantities or assign conflicting values to the same parameter, one definition may overwrite the other during configuration or model initialization.

This can lead to unintended parameter values and **undefined or unpredictable simulation behaviour**.

In particular:

- Ensure that models sharing parameters use compatible definitions and configurations.
- Verify that all related models use the intended configuration.

For example:
`BYCharm` must use the same configuration as `KNOTunedQPMDISPXSec` as they share common variables. If different configurations are used, for example `KNOTunedQPMDISPXSec/Default` in combination with `BYCharm/BY21`, the resulting behavoir is undefined.


Whenever a model is added or modified, its parameter names and configuration dependencies must be checked against the other models in the CMC.

**The key requirement is that every parameter name must refer to one well-defined physical quantity and have a consistent configuration throughout the tune.**

### Comprehensive Model Configuration

A CMC is identified by:

`Gdd_MMv`

| Field | Description                                       |
| ----- | ------------------------------------------------- |
| `G`   | Tune author (`G` denotes the GENIE collaboration) |
| `dd`  | Year in which the configuration was developed     |
| `MM`  | Model configuration branch identifier             |
| `v`   | Variation or branch version (`a`, `b`, `c`, ...)  |

Test configurations use the prefix `GTEST` instead of `G`. These are temporary configurations and are not official GENIE tunes.

### Tune

A tune is identified by:

`Gdd_MMv_PP_xxx`

| Field     | Description                                                          |
| --------- | -------------------------------------------------------------------- |
| `Gdd_MMv` | Underlying comprehensive model configuration                         |
| `PP`      | Tuned parameter set identifier                                       |
| `xxx`     | Tuning dataset identifier, including the associated weighting scheme |

The identifiers `00` for the parameter set and dataset indicate that the configuration has not been tuned by GENIE.
