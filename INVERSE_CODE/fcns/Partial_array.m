classdef Partial_array < handle
    properties
        size;
        len;
        gpuData;
        maplist;
        slice0;
        func;
        % ext_param;
        index;

    end
    
    methods
        function [this, final_slice] = Partial_array(xy_size, len, datatype, func, slice0)
            size = zeros(1, 3, "uint16");
            size(1:2) = xy_size;
            size(3) = 2+floor(sqrt(len*2));
            this.size = size;
            this.len = len;
            if datatype == "complex"
                this.gpuData = complex(zeros(size,"gpuArray"));
            else
                this.gpuData = zeros(size, datatype,"gpuArray");
            end
            this.func = func;
            this.slice0 = slice0;
            this.maplist = zeros(size(1), 1, "uint16");
            final_slice = this.compute();
        end
        
        function fi = compute(this)
            this.index = 0;
            max_skip = floor(sqrt(this.len*2));
            skip = 0;
            curr = this.slice0;
            
            for i=1:this.len
                curr = this.func(curr, i);
                if skip==0
                    this.index = this.index+1;
                    if this.index>this.size(3)
                        error("internal error: preallocated size too small");
                    end
                    
                    this.gpuData(:,:,this.index) = curr;
                    this.maplist(this.index) = i;
                    skip = max_skip;
                    max_skip = max(max_skip-1, 0);
                else
                    skip = skip - 1;
                end
            end
            this.index = this.index+1;
            this.maplist(this.index) = Inf;
            fi = curr;
        end
        
        function data = get(this, i)
            if i > this.maplist(this.index)
                error("get index cannot increase")
            end
            % find the closest recorded slice *before* i
            while i < this.maplist(this.index)
                this.index = this.index - 1;
            end
            % re-calculate until i
            while i > this.maplist(this.index)
                this.index = this.index + 1;
                if this.index>this.size(3)
                    error("internal error: preallocated size too small");
                end
                this.maplist(this.index) = this.maplist(this.index-1) + 1;
                this.gpuData(:,:,this.index) = this.func(...
                    this.gpuData(:,:,this.index-1), ...
                    this.maplist(this.index));
            end
            
            data = this.gpuData(:,:,this.index);
        end
            
    end
    
end