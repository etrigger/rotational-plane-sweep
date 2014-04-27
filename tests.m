clear; close all; clc;

%vertices = environment();

vertices = load('vertices_book.mat');
vertices = vertices.vertices

[ edges ] = RPS( vertices )