// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_EXTRASBUNDLE_HPP
#define ZENTAS_EXTRASBUNDLE_HPP

#include "claransextrasbundle.hpp"


class ExtrasBundle{

public:
  ClaransExtrasBundle clarans;
  ExtrasBundle(size_t max_proposals, size_t patient):clarans(max_proposals, patient){}
  
};


#endif
