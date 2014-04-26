clear; close all; clc;

%vertices = environment();

vertices = load('vertices2.mat');
vertices = vertices.vertices

[ edges ] = RPS( vertices )