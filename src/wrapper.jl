struct LLGWrapper <: OrdinaryLLGAlgorithm end

isfsal(alg::LLGEulerHeun) = false
alg_order(alg::LLGEulerHeun) = 1
