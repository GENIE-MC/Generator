#include "Physics/QuasiElastic/XSection/ManualResponseTensor.h"

genie::ManualResponseTensor::ManualResponseTensor(const std::complex<double> in_RespArray[4][4]
	)
{
	memcpy(RespTensor, in_RespArray, sizeof(RespTensor)); //Copies input array and sets e
}

std::complex<double> genie::ManualResponseTensor::operator()(TensorIndex_t mu, TensorIndex_t nu)
  	const /*override*/
{
	std::complex<double> result = this->RespTensor[mu][nu];

	return result;
}
