# Network Formulations

## Type Hierarchy
We begin with the top of the hierarchy, where we can distinguish between conic and non-conic power flow models.
```julia
AbstractConicForms = Union{AbstractConicPowerFormulation, AbstractBFConicForm}
AbstractConicPowerFormulation <: PowerModels.AbstractPowerFormulation
AbstractBFForm <: PowerModels.AbstractPowerFormulation
AbstractBFQPForm <: PowerModels.AbstractBFForm
AbstractBFConicForm <: PowerModels.AbstractBFForm
```

We begin with the top of the hierarchy, where we can distinguish between AC and DC power flow models.
```julia
AbstractACPForm <: PowerModels.AbstractPowerFormulation
AbstractDCPForm <: PowerModels.AbstractPowerFormulation
AbstractNLPUBFForm <: PowerModels.AbstractBFQPForm
AbstractConicUBFForm <: PowerModels.AbstractBFConicForm
AbstractLPUBFForm <: AbstractNLPUBFForm
```

From there, different forms are possible:
```julia
StandardACPForm <: PowerModels.AbstractACPForm
StandardDCPForm <: PowerModels.AbstractDCPForm

SDPUBFForm <: AbstractConicUBFForm
SOCNLPUBFForm <: AbstractNLPUBFForm
SOCConicUBFForm <: AbstractConicUBFForm

LPfullUBFForm <: AbstractLPUBFForm
LPdiagUBFForm <: AbstractLPUBFForm
```



## Power Models
Each of these forms can be used as the type parameter for a PowerModel:
```julia
ACPPowerModel = GenericPowerModel{PowerModels.StandardACPForm}
DCPPowerModel = GenericPowerModel{PowerModels.StandardDCPForm}

SDPUBFPowerModel = GenericPowerModel{SDPUBFForm}
SOCNLPUBFPowerModel = GenericPowerModel{SOCNLPUBFForm}
SOCConicUBFPowerModel = GenericPowerModel{SOCConicUBFForm}

LPfullUBFPowerModel = GenericPowerModel{LPfullUBFForm}
LPdiagUBFPowerModel = GenericPowerModel{LPdiagUBFForm}
LPLinUBFPowerModel = PMs.GenericPowerModel{LPLinUBFForm}
```


## Union Types

To support both conic and quadratically-constrained formulation variants for the unbalanced branch flow model, the union type `AbstractUBFForm` is defined. These formulations extend `AbstractBFForm` and are therefore also `AbstractWForms`.

```julia
AbstractUBFForm = Union{AbstractNLPUBFForm, AbstractConicUBFForm}
AbstractWForms = Union{AbstractWRForms, AbstractBFForm}
```