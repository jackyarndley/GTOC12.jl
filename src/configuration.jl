abstract type CoordinateConfig end

struct CartesianConfig <: CoordinateConfig end
struct ModifiedEquinoctialConfig <: CoordinateConfig end
struct ModifiedOrbitalConfig <: CoordinateConfig end

abstract type ThrustConfig end

struct ImpulsiveConfig <: ThrustConfig end
struct ZeroOrderHoldConfig <: ThrustConfig end
struct FirstOrderHoldConfig <: ThrustConfig end

abstract type ObjectiveConfig end

struct VelocityConfig <: ObjectiveConfig end
struct LoggedMassConfig <: ObjectiveConfig end
struct FinalConfig <: ObjectiveConfig end
struct MassConfig <: ObjectiveConfig end