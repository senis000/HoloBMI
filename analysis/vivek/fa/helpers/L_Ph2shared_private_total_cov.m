function [shared_cov, private_cov, total_cov] = ...
    L_Ph2shared_private_total_cov(L, Ph)
%4.10.17
%takes in parameters from 'fastfa.m', returns the shared, private and total
%covariance matrices
%L - FA factors
%Ph - vector, diagonal of the private covariance matrix

shared_cov      = L*L'; 
private_cov     = diag(Ph); 
total_cov       = shared_cov + private_cov;