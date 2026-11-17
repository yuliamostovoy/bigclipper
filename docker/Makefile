VERSION = 2.2
TAG1 = quay.io/ymostovoy/clipper:$(VERSION)
TAG2 = quay.io/ymostovoy/clipper:latest

all: build push

# build context is the repository root so the image gets the local scripts
build:
	docker build -f Dockerfile -t $(TAG1) -t $(TAG2) ..

push:
	docker push $(TAG1)
	docker push $(TAG2)
