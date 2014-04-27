clear; close all; clc;

%vertices = environment();

vertices = load('vertices_lab.mat');
vertices = vertices.vertices

[ edges ] = RPS( vertices )