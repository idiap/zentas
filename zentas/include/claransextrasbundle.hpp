// Copyright (c) 2016 Idiap Research Institute, http://www.idiap.ch/
// Written by James Newling <jnewling@idiap.ch>

#ifndef ZENTAS_CLARANSEXTRASBUNDLE_HPP
#define ZENTAS_CLARANSEXTRASBUNDLE_HPP

class ClaransExtrasBundle{

public:
  size_t max_proposals;
  bool patient;
  
  ClaransExtrasBundle(size_t max_proposals_, bool patient_):max_proposals(max_proposals_), patient(patient_) {}
};


#endif
