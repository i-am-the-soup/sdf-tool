const std = @import("std");
const stb = @cImport({
    @cInclude("stb_image.h");
    @cInclude("stb_image_write.h");
});

const Point = struct {
    dx: i32,
    dy: i32,

    fn min(a: Point, b: Point) Point {
        if (a.distance2() < b.distance2()) {
            return a;
        } else {
            return b;
        }
    }

    fn distance2(self: Point) i32 {
        return self.dx *| self.dx +| self.dy *| self.dy;
    }
};

const Grid = struct {
    width: usize,
    height: usize,
    grid: []Point,

    fn init(image: []const u8, width: usize, height: usize, invert: bool, grid: []Point) Grid {
        std.debug.assert(image.len == width * height);

        for (image, grid) |i, *o| {
            if (!invert) {
                if (i > std.math.maxInt(u8) / 2) {
                    o.* = .{ .dx = 0, .dy = 0 };
                } else {
                    o.* = .{ .dx = std.math.maxInt(i16), .dy = std.math.maxInt(i16) };
                }
            } else {
                if (i > std.math.maxInt(u8) / 2) {
                    o.* = .{ .dx = std.math.maxInt(i16), .dy = std.math.maxInt(i16) };
                } else {
                    o.* = .{ .dx = 0, .dy = 0 };
                }
            }
        }

        return Grid{
            .width = width,
            .height = height,
            .grid = grid,
        };
    }

    fn process(self: *Grid) void {
        for (0..self.height) |y| {
            for (0..self.width) |x| {
                var p = self.get(x, y);
                p = Point.min(p, self.getOffset(x -% 1, y -% 1, -1, -1));
                p = Point.min(p, self.getOffset(x, y -% 1, 0, -1));
                p = Point.min(p, self.getOffset(x +% 1, y -% 1, 1, -1));
                p = Point.min(p, self.getOffset(x -% 1, y, -1, 0));
                self.put(x, y, p);
            }

            for (0..self.width) |x0| {
                const x = (self.width - 1) - x0;
                var p = self.get(x, y);
                p = Point.min(p, self.getOffset(x +% 1, y, 1, 0));
                self.put(x, y, p);
            }
        }

        for (0..self.height) |y0| {
            const y = (self.height - 1) - y0;
            for (0..self.width) |x| {
                var p = self.get(x, y);
                p = Point.min(p, self.getOffset(x +% 1, y +% 1, 1, 1));
                p = Point.min(p, self.getOffset(x, y +% 1, 0, 1));
                p = Point.min(p, self.getOffset(x -% 1, y +% 1, -1, 1));
                p = Point.min(p, self.getOffset(x +% 1, y, 1, 0));
                self.put(x, y, p);
            }

            for (0..self.width) |x0| {
                const x = (self.width - 1) - x0;
                var p = self.get(x, y);
                p = Point.min(p, self.getOffset(x -% 1, y, -1, 0));
                self.put(x, y, p);
            }
        }
    }

    fn writeSdf(
        self: Grid,
        out_buf: []f32,
    ) void {
        std.debug.assert(out_buf.len == self.width * self.height);

        for (self.grid, out_buf) |i, *o| {
            const d = @sqrt(@intToFloat(f32, i.distance2())) / 2048.0 / @sqrt(2.0);
            o.* = d;
        }
    }

    fn get(self: Grid, x: usize, y: usize) Point {
        if (x >= self.width or y >= self.height) {
            return Point{ .dx = std.math.maxInt(i16), .dy = std.math.maxInt(i16) };
        }

        return self.grid[y * self.width + x];
    }

    fn getOffset(self: Grid, x: usize, y: usize, dx: i32, dy: i32) Point {
        var p = self.get(x, y);
        p.dx +|= dx;
        p.dy +|= dy;
        return p;
    }

    fn put(self: *Grid, x: usize, y: usize, p: Point) void {
        self.grid[y * self.width + x] = p;
    }
};

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    defer if (gpa.deinit())
        std.log.err("Memory leaked on exit!", .{});
    const allocator = gpa.allocator();

    var maybe_input_dir_path: ?[]const u8 = null;
    var maybe_output_file_path: ?[]const u8 = null;
    var maybe_blending_steps: ?usize = 20;

    var args_it = try std.process.argsWithAllocator(allocator);
    defer args_it.deinit();

    _ = args_it.next() orelse unreachable;

    while (args_it.next()) |arg| {
        if (maybe_input_dir_path == null) {
            maybe_input_dir_path = arg;
        } else if (maybe_output_file_path == null) {
            maybe_output_file_path = arg;
        } else {
            maybe_blending_steps = std.fmt.parseInt(usize, arg, 10) catch null;
        }
    }

    if (maybe_input_dir_path == null or maybe_output_file_path == null or maybe_blending_steps == null) {
        std.log.err("Invalid usage!\nUsage: sdf_tool <INPUT_FILE_DIRECTORY> <OUTPUT_FILE_NAME> [BLENDING_STEPS (default 20)]", .{});
        return error.InvalidUsage;
    }

    const input_dir_path = maybe_input_dir_path.?;
    const output_file_path = try std.cstr.addNullByte(allocator, maybe_output_file_path.?);
    defer allocator.free(output_file_path);

    var input_file_paths = try findInputFilePaths(allocator, input_dir_path);
    defer {
        for (input_file_paths) |path| {
            allocator.free(path);
        }
        allocator.free(input_file_paths);
    }

    if (input_file_paths.len < 2) {
        std.log.err("Expected two or more input images!", .{});
        return error.NotEnoughImages;
    }

    var width: usize = 0;
    var height: usize = 0;

    try validateInputFiles(input_file_paths, &width, &height);

    var images = try allocator.alloc([]u8, input_file_paths.len);
    defer allocator.destroy(images);

    for (input_file_paths, images) |path, *img| {
        img.* = try allocator.alloc(u8, width * height);
        try loadImage(path, img.*);
    }
    defer for (images) |img| allocator.free(img);

    var output_data = try allocator.alloc(f32, width * height);
    defer allocator.free(output_data);
    std.mem.set(f32, output_data, 0.0);

    var grid_data = try allocator.alloc(Point, width * height);
    defer allocator.free(grid_data);

    var sdf_buffer = try allocator.alloc(f32, width * height * 2);
    defer allocator.free(sdf_buffer);

    var sdf1_outer = try allocator.alloc(f32, width * height);
    defer allocator.free(sdf1_outer);

    var sdf1_inner = try allocator.alloc(f32, width * height);
    defer allocator.free(sdf1_inner);

    var sdf2_outer = try allocator.alloc(f32, width * height);
    defer allocator.free(sdf2_outer);

    var sdf2_inner = try allocator.alloc(f32, width * height);
    defer allocator.free(sdf2_inner);

    for (images[0 .. images.len - 1], 0..) |_, i| {
        {
            var grid_outer = Grid.init(images[i], width, height, false, grid_data);
            grid_outer.process();
            grid_outer.writeSdf(sdf1_outer);

            var grid_inner = Grid.init(images[i], width, height, true, grid_data);
            grid_inner.process();
            grid_inner.writeSdf(sdf1_inner);
        }
        {
            var grid_outer = Grid.init(images[i + 1], width, height, false, grid_data);
            grid_outer.process();
            grid_outer.writeSdf(sdf2_outer);

            var grid_inner = Grid.init(images[i + 1], width, height, true, grid_data);
            grid_inner.process();
            grid_inner.writeSdf(sdf2_inner);
        }

        const blending_steps = maybe_blending_steps.?;
        for (0..blending_steps) |j| {
            const t = @intToFloat(f32, j) / @intToFloat(f32, blending_steps);
            for (sdf1_outer, sdf1_inner, sdf2_outer, sdf2_inner, output_data) |a_out, a_in, b_out, b_in, *o| {
                const s = (1.0 - t) * (a_out - a_in) + t * (b_out - b_in);
                if (s <= 0) {
                    const addend = 1.0 / @intToFloat(f32, blending_steps * (images.len - 1));
                    o.* += addend;
                }
            }
        }
    }

    // normalize image
    {
        var max: f32 = 0.0;
        for (output_data) |o| {
            if (o > max) {
                max = o;
            }
        }

        for (output_data) |*o| {
            o.* /= max;
        }
    }

    var output_image = try allocator.alloc(u8, width * height);
    defer allocator.free(output_image);

    for (output_data, output_image) |i, *o| {
        o.* = @truncate(u8, @min(255, @floatToInt(usize, 255.0 * i)));
    }

    if (stb.stbi_write_png(
        output_file_path.ptr,
        @intCast(c_int, width),
        @intCast(c_int, height),
        1,
        output_image.ptr,
        @intCast(c_int, width),
    ) == 0) {
        std.log.err("Failed to write output file '{s}'", .{output_file_path});
        return error.OutputWriteFailed;
    }
}

fn findInputFilePaths(allocator: std.mem.Allocator, dir_path: []const u8) ![][:0]const u8 {
    var dir = try std.fs.cwd().openIterableDir(dir_path, .{});
    defer dir.close();

    var input_file_paths = std.ArrayList([:0]const u8).init(allocator);
    defer input_file_paths.deinit();

    var it = dir.iterate();
    while (try it.next()) |file| {
        if (file.kind != .File) {
            continue;
        }

        var path = try dir.dir.realpathAlloc(allocator, file.name);
        defer allocator.free(path);

        try input_file_paths.append(try std.cstr.addNullByte(allocator, path));
    }

    std.sort.sort([]const u8, input_file_paths.items, {}, stringAsc);

    return try input_file_paths.toOwnedSlice();
}

fn stringAsc(context: void, a: []const u8, b: []const u8) bool {
    _ = context;
    const min_len = @min(a.len, b.len);
    for (a[0..min_len], b[0..min_len]) |x, y| {
        if (x == y) {
            continue;
        } else {
            return x < y;
        }
    }

    return false;
}

fn validateInputFiles(input_file_paths: [][:0]const u8, out_width: *usize, out_height: *usize) !void {
    for (input_file_paths) |path| {
        var w: c_int = undefined;
        var h: c_int = undefined;
        var c: c_int = undefined;

        if (stb.stbi_info(path.ptr, &w, &h, &c) == 0) {
            std.log.err("Couldn't fetch image info for '{s}'! Make sure the input folder consits of only images files!", .{path});
            return error.ImageInfoFetchFailed;
        }

        if (c != 1) {
            std.log.err("Expected single-channel image! Please convert your images to a grayscale color space!", .{});
            return error.NonGrayscaleImage;
        }

        if (out_width.* == 0 or out_height.* == 0) {
            out_width.* = @intCast(usize, w);
            out_height.* = @intCast(usize, h);
        }

        if (w != out_width.* or h != out_height.*) {
            std.log.err("Inconsitent image dimensions! Make sure all images are same the same dimensions!", .{});
            return error.InconsistentImageDimensions;
        }
    }
}

fn loadImage(path: [:0]const u8, out_buf: []u8) !void {
    var w: c_int = undefined;
    var h: c_int = undefined;
    var c: c_int = undefined;

    var data_p = stb.stbi_load(path.ptr, &w, &h, &c, 0);
    if (data_p == null) {
        std.log.err("Failed to load image: '{s}'", .{path});
        return error.ImageLoadError;
    }
    defer std.c.free(data_p);

    // All of these should hold if validteInputFiles passes.
    std.debug.assert(c == 1);
    std.debug.assert(out_buf.len == w * h);

    std.mem.copy(u8, out_buf, data_p[0..@intCast(usize, w * h)]);
}
